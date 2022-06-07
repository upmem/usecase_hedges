#ifndef __XFERITF__
#define __XFERITF__

#ifndef DPU_ENV
#include <dpu>
#include <cassert>
#endif

#include <common.h>

enum SharedObjType {
    UcharType,
    CharType,
    UllongType,
    IntType,
    VecCharType,
    VecUcharType,
    VecIntType,
    VecUllongType,
};

enum DpuAlloctype { none, mem, buddy };

#ifndef DPU_ENV

namespace dpu {
template <typename T> class Tensor2d {
public:
    Tensor2d(uint64_t nr_dpus, uint64_t shapes_[])
        : nr_dpus(nr_dpus)
    {
        rank = 2;
        shapes[0] = shapes_[0];
        shapes[1] = shapes_[1];
        shapes_aligned[0] = shapes_[0];
        shapes_aligned[1] = XFER_MEM_ALIGN(shapes_[1] * sizeof(T)) / sizeof(T);
        aligned_size = shapes_aligned[0] * shapes_aligned[1];

        data = new std::vector<std::vector<T>>(nr_dpus, std::vector<T>(aligned_size));
    }

    T &operator()(uint64_t dpu_index, uint64_t first, uint64_t second)
    {
        assert(dpu_index < nr_dpus);
        return (*data)[dpu_index][first * shapes_aligned[1] + second];
    }

    ~Tensor2d() { delete data; }

    uint64_t nr_dpus;
    uint64_t aligned_size;
    uint8_t rank;
    uint64_t shapes[2];
    uint64_t shapes_aligned[2];
    std::vector<std::vector<T>> *data;
};

/*
 * generic class for HOST/DPU data exchange with arbitrary type
 * */
class xferItf {
public:
    using DescriptorImpl = uint64_t;

    static uint64_t DpuAllign(uint64_t val) { return XFER_MEM_ALIGN(val); }

    xferItf(const char *buffSymbolName, dpu_set_t &set)
        : buffOffset(0)
        , buffSymbolName(buffSymbolName)
        , dpu_set(&set)
        , async(false)
    {
        DPU_ASSERT(dpu_get_nr_dpus(set, &nr_dpus));
        flags = async ? DPU_XFER_ASYNC : DPU_XFER_DEFAULT;
    }

    void offsetInc(DescriptorImpl &offset) { buffOffset += offset; }

    void save() { savedBuffOffset = buffOffset; }

    void restore() { buffOffset = savedBuffOffset; }

    DescriptorImpl offset() { return buffOffset; }

    void broadcastDesc(DescriptorImpl type)
    {
        auto sz = XFER_MEM_ALIGN(sizeof(DescriptorImpl));
        DPU_ASSERT(dpu_broadcast_to(*dpu_set, buffSymbolName, offset(), (uint8_t *)(&type), sz, flags));
        offsetInc(sz);
    }

    void pushDesc(DescriptorImpl desc)
    {
        auto sz = XFER_MEM_ALIGN(sizeof(DescriptorImpl));
        DPU_FOREACH (*dpu_set, dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, (uint8_t *)(&desc)));
        }
        DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_TO_DPU, buffSymbolName, offset(), sz, flags));
        offsetInc(sz);
    }

    void getDesc(std::vector<DescriptorImpl> &desc)
    {
        uint32_t each_dpu;
        auto sz = XFER_MEM_ALIGN(sizeof(DescriptorImpl));
        DPU_FOREACH (*dpu_set, dpu, each_dpu) {
            DPU_ASSERT(dpu_prepare_xfer(dpu, (uint8_t *)(&(desc.data()[each_dpu]))));
        }
        DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_FROM_DPU, buffSymbolName, offset(), sz, flags));
        offsetInc(sz);
    }

    template <typename T> void scalarBroadcast(SharedObjType objType, T val)
    {
        broadcastDesc(objType);
        broadcastDesc(1);

        DescriptorImpl packetSz = XFER_MEM_ALIGN(sizeof(T));
        DPU_ASSERT(dpu_broadcast_to(*dpu_set, buffSymbolName, offset(), (uint8_t *)(&val), packetSz, flags));
        offsetInc(packetSz);
    }

    template <typename T> void vectorBroadcast(SharedObjType objType, std::vector<T> &val)
    {
        auto len = val.size();
        auto *data = val.data();

        broadcastDesc(objType);
        broadcastDesc(len);

        DescriptorImpl packetSz = XFER_MEM_ALIGN(len * sizeof(T));
        DPU_ASSERT(dpu_broadcast_to(*dpu_set, buffSymbolName, offset(), (uint8_t *)(data), packetSz, DPU_XFER_DEFAULT));
        offsetInc(packetSz);
    }

    template <typename T> void vectorBroadcast(SharedObjType objType, NRvector<T> &val)
    {
        auto len = val.size();
        auto *data = val.data();
        pushDesc(objType);
        pushDesc(len);
        DescriptorImpl packetSz = XFER_MEM_ALIGN(len * sizeof(T));
        DPU_ASSERT(dpu_broadcast_to(*dpu_set, buffSymbolName, offset(), (uint8_t *)(data), packetSz, DPU_XFER_DEFAULT));
        offsetInc(packetSz);
    }

    template <typename T> void vectorPush(SharedObjType objType, std::vector<std::vector<T>> &val)
    {
        uint32_t each_dpu;

        auto dpu_dim = val.size();
        assert((dpu_dim == nr_dpus) && "incompatible destination vector capacity");
        auto data_dim = val[0].size();
        for (auto &data : val)
            assert(data.size() == data_dim && "vector for DPUs must have equal size");

        DescriptorImpl packetSz = XFER_MEM_ALIGN(data_dim * sizeof(T));

        // NOTE objtypes are not checked yet
        pushDesc(objType);
        pushDesc(data_dim);

        DPU_FOREACH (*dpu_set, dpu, each_dpu) {
            auto *dpuData = val[each_dpu].data();
            DPU_ASSERT(dpu_prepare_xfer(dpu, dpuData));
        }
        DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_TO_DPU, buffSymbolName, offset(), packetSz, flags));
        offsetInc(packetSz);
    }

    template <typename T> void vectorGet(std::vector<std::vector<T>> &val)
    {
        uint32_t each_dpu;
        auto dpu_dim = val.capacity();
        assert((dpu_dim == nr_dpus) && "incompatible destination vector capacity");
        auto data_dim = val[0].capacity();

        std::vector<DescriptorImpl> rxObjTypes(nr_dpus);
        std::vector<DescriptorImpl> rxObjLens(nr_dpus);

        // NOTE objtypes are not checked yet
        getDesc(rxObjTypes);
        getDesc(rxObjLens);

        for (auto &rxlen : rxObjLens)
            assert((rxlen == data_dim) && "incompatible destination vector capacity");

        DescriptorImpl packetSz = XFER_MEM_ALIGN(data_dim * sizeof(T));

        DPU_FOREACH (*dpu_set, dpu, each_dpu) {
            auto *dpuData = val[each_dpu].data();
            DPU_ASSERT(dpu_prepare_xfer(dpu, dpuData));
        }
        DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_FROM_DPU, buffSymbolName, offset(), packetSz, flags));
        offsetInc(packetSz);
    }

    template <class T> void broadcast(T) = delete;
    template <class T, class U> void broadcast(T, U) = delete;
    template <class T> void get(T) = delete;
    template <class T, class U> void get(T, U) = delete;

    void broadcast(Ullong val) { scalarBroadcast(UllongType, val); }

    void broadcast(Uchar val) { scalarBroadcast(UcharType, val); }

    void broadcast(Int val) { scalarBroadcast(IntType, val); }

    template <typename T> void push(Tensor2d<T> &t)
    {
        uint32_t each_dpu;

        assert((t.nr_dpus == nr_dpus) && "incompatible destination vector capacity");

        DescriptorImpl packetSz = t.aligned_size * sizeof(T);

        // NOTE objtypes are not checked yet
        pushDesc(CharType);
        pushDesc(t.shapes[0]);
        pushDesc(t.shapes[1]);

        DPU_FOREACH (*dpu_set, dpu, each_dpu) {
            auto *dpuData = (*t.data)[each_dpu].data();
            DPU_ASSERT(dpu_prepare_xfer(dpu, dpuData));
        }
        DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_TO_DPU, buffSymbolName, offset(), packetSz, flags));
        offsetInc(packetSz);
    }

    template <typename T> void get(Tensor2d<T> &t)
    {
        uint32_t each_dpu;
        auto dpu_dim = (*t.data).size();
        assert((dpu_dim == nr_dpus) && "incompatible destination vector capacity");

        std::vector<DescriptorImpl> rxObjTypes(nr_dpus);
        std::vector<DescriptorImpl> shapes_0(nr_dpus);
        std::vector<DescriptorImpl> shapes_1(nr_dpus);
        // NOTE objtypes are not checked yet
        getDesc(rxObjTypes);
        getDesc(shapes_0);
        getDesc(shapes_1);
        for (auto &shape : shapes_0) {
            assert((shape == t.shapes[0]) && "incompatible destination vector capacity");
        }
        for (auto &shape : shapes_1)
            assert((shape == t.shapes[1]) && "incompatible destination vector capacity");

        DescriptorImpl packetSz = t.aligned_size * sizeof(T);

        DPU_FOREACH (*dpu_set, dpu, each_dpu) {
            auto *dpuData = (*t.data)[each_dpu].data();
            DPU_ASSERT(dpu_prepare_xfer(dpu, dpuData));
        }
        DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_FROM_DPU, buffSymbolName, offset(), packetSz, flags));
        offsetInc(packetSz);
    }

    void push(std::vector<std::vector<char>> &val) { vectorPush(CharType, val); }

    void push(std::vector<std::vector<Uchar>> &val) { vectorPush(UcharType, val); }

    void push(std::vector<std::vector<Ullong>> &val) { vectorPush(UllongType, val); }

    void push(std::vector<std::vector<Int>> &val) { vectorPush(IntType, val); }

    void broadcast(VecUchar val) { vectorBroadcast(VecUcharType, val); }

    void broadcast(VecInt val) { vectorBroadcast(VecIntType, val); }

    void broadcast(std::vector<std::vector<Int>> &val) { vectorBroadcast(IntType, val); }

    void broadcast(VecUllong val) { vectorBroadcast(VecUllongType, val); }

    void get(std::vector<std::vector<Uchar>> &val) { vectorGet(val); }

    void get(std::vector<std::vector<Ullong>> &val) { vectorGet(val); }

private:
    size_t buffOffset;
    size_t savedBuffOffset;
    const char *buffSymbolName;
    struct dpu_set_t *dpu_set;
    dpu_set_t dpu;
    uint32_t nr_dpus;
    bool async;
    dpu_xfer_flags_t flags;
};
} // end namespace dpu
#else

typedef uint64_t DescriptorImpl;

typedef struct Tensor2d {
    __mram_ptr uint8_t **mram_addr;
    uint64_t shapes[2];
    uint64_t aligned_shapes[2];
    DescriptorImpl dataType;

#undef DIM_MAX
} Tensor2d;

void check2dTensorUnsigned(Tensor2d *t)
{
    printf("[check2dInputTensor]\n");
    printf("  shape[0] %lu \n", t->shapes[0]);
    printf("  shape[1] %lu \n", t->shapes[1]);
    printf("  align shape[0] %lu \n", t->aligned_shapes[0]);
    printf("  align shape[1] %lu \n", t->aligned_shapes[1]);
    uint64_t wram_size = sizeof(uint32_t) * t->aligned_shapes[1];
    __dma_aligned uint32_t *loc = buddy_alloc(wram_size);
    for (uint64_t i = 0; i < t->shapes[0]; i++) {
        mram_read(t->mram_addr[i], loc, wram_size);
        for (uint64_t j = 0; j < t->shapes[1]; j++)
            printf("  data[%lu][%lu] %u \n", i, j, (unsigned)(loc[j]));
    }
    buddy_free(loc);
}
void check2dTensorUchar(Tensor2d *t)
{
    printf("[check2dInputTensor]\n");
    printf("  shape[0] %lu \n", t->shapes[0]);
    printf("  shape[1] %lu \n", t->shapes[1]);
    printf("  align shape[0] %lu \n", t->aligned_shapes[0]);
    printf("  align shape[1] %lu \n", t->aligned_shapes[1]);
    uint64_t wram_size = sizeof(Uchar) * t->aligned_shapes[1];
    printf("  wram size %lu \n", wram_size);
    Uchar loc[400];
    for (uint64_t i = 0; i < t->shapes[0]; i++) {
        mram_read(t->mram_addr[i], loc, wram_size);
        for (uint64_t j = 0; j < t->shapes[1]; j++)
            printf("  data[%lu][%lu] %u \n", i, j, (unsigned)(loc[j]));
    }
    printf("  wram size 1 %lu \n", wram_size);
}

void free2dTensor(Tensor2d *t) { buddy_free(t->mram_addr); }

typedef struct xferItf {
    __mram_ptr uint8_t *buffer;
    __dma_aligned DescriptorImpl offset, savedOffset;

    void (*init)(struct xferItf *, __mram_ptr uint8_t *);
    void (*save)(struct xferItf *, uint64_t *);
    void (*restore)(struct xferItf *, uint64_t *);
    void (*offsetInc)(struct xferItf *, DescriptorImpl);
    void (*getDesc)(struct xferItf *, DescriptorImpl *);
    void (*pushDesc)(struct xferItf *, DescriptorImpl);
    void (*getUchar)(struct xferItf *, Uchar *);
    void (*getUllong)(struct xferItf *, Ullong *);
    void (*getInt)(struct xferItf *, Int *);
    void (*getVecUchar)(struct xferItf *, VecUchar *, DescriptorImpl *, enum DpuAlloctype);
    void (*getVecChar)(struct xferItf *, char **, DescriptorImpl *, enum DpuAlloctype);
    void (*getVecInt)(struct xferItf *, VecInt *, DescriptorImpl *, enum DpuAlloctype);
    void (*getVecUllong)(struct xferItf *, VecUllong *, DescriptorImpl *, enum DpuAlloctype);
    void (*pushVecUchar)(struct xferItf *, VecUchar, DescriptorImpl);
    void (*getVecCharNoAlloc)(struct xferItf *, DescriptorImpl *);
    void (*getTensor2dUnsigned)(struct xferItf *, Tensor2d *);
    void (*pushTensor2dUnsigned)(struct xferItf *, Tensor2d *, uint64_t[]);
    void (*getTensor2dUchar)(struct xferItf *, Tensor2d *);
    void (*pushTensor2dUchar)(struct xferItf *, Tensor2d *, uint64_t[]);
    void (*pushTensor2dUINT64)(struct xferItf *, Tensor2d *, uint64_t[]);
} xferItf;

#define XFERITF_INIT()                                                                                                           \
    {                                                                                                                            \
        .init = xferItf_init, .save = xferItf_save, .restore = xferItf_restore, .offsetInc = xferItf_offsetInc,                  \
        .getDesc = xferItf_getDesc, .pushDesc = xferItf_pushDesc, .getUchar = xferItf_getUchar, .getUllong = xferItf_getUllong,  \
        .getInt = xferItf_getInt, .getVecUchar = xferItf_getVecUchar, .getVecChar = xferItf_getVecChar,                          \
        .getVecInt = xferItf_getVecInt, .getVecUllong = xferItf_getVecUllong, .pushVecUchar = xferItf_pushVecUchar,              \
        .getTensor2dUnsigned = xferItf_getTensor2dUnsigned, .pushTensor2dUnsigned = xferItf_pushTensor2dUnsigned,                \
        .getTensor2dUchar = xferItf_getTensor2dUchar, .pushTensor2dUchar = xferItf_pushTensor2dUchar,                            \
        .pushTensor2dUINT64 = xferItf_pushTensor2dUINT64,                                                                        \
    }

#define VECTOR_ALLOCATE_AND_FETCH_MRAM_TO_WRAM(VECVALPTR, VECLEN, TYPE, ATYPE)                                                   \
    __dma_aligned DescriptorImpl type;                                                                                           \
    this->getDesc(this, &type);                                                                                                  \
    this->getDesc(this, VECLEN);                                                                                                 \
    DescriptorImpl sz = XFER_MEM_ALIGN(*VECLEN * sizeof(TYPE));                                                                  \
    if (ATYPE == buddy) {                                                                                                        \
        *VECVALPTR = buddy_alloc(sz);                                                                                            \
    } else if (ATYPE == mem) {                                                                                                   \
        *VECVALPTR = mem_alloc(sz);                                                                                              \
    }                                                                                                                            \
    mram_read(this->buffer + this->offset, *VECVALPTR, sz);                                                                      \
    this->offsetInc(this, sz);

#define VECTOR_FETCH_WRAM_TO_MRAM(VECVALPTR, VECLEN, TYPE, TYPEID)                                                               \
    this->pushDesc(this, TYPEID);                                                                                                \
    this->pushDesc(this, VECLEN);                                                                                                \
    DescriptorImpl sz = XFER_MEM_ALIGN(VECLEN * sizeof(TYPE));                                                                   \
    mram_write(VECVALPTR, this->buffer + this->offset, sz);                                                                      \
    this->offsetInc(this, sz);

#define SCALAR_FETCH_MRAM_TO_WRAM(VALPTR, TYPE)                                                                                  \
    __dma_aligned DescriptorImpl type, len;                                                                                      \
    this->getDesc(this, &type);                                                                                                  \
    this->getDesc(this, &len);                                                                                                   \
    DescriptorImpl sz = XFER_MEM_ALIGN(len * sizeof(TYPE));                                                                      \
    mram_read(this->buffer + this->offset, VALPTR, sz);                                                                          \
    this->offsetInc(this, sz);

void xferItf_save(xferItf *this, uint64_t *reg) { *reg = this->offset; }
void xferItf_restore(xferItf *this, uint64_t *reg) { this->offset = *reg; }

void xferItf_init(xferItf *this, __mram_ptr uint8_t *buffer)
{
    this->buffer = buffer;
    this->offset = 0;
}

void xferItf_offsetInc(xferItf *this, DescriptorImpl offset) { this->offset += offset; }

void xferItf_getDesc(xferItf *this, DescriptorImpl *desc)
{
    __dma_aligned DescriptorImpl sz = sizeof(DescriptorImpl);
    mram_read(this->buffer + this->offset, desc, sz);
    this->offsetInc(this, sz);
}

void xferItf_pushDesc(xferItf *this, DescriptorImpl desc)
{
    __dma_aligned DescriptorImpl sz = sizeof(DescriptorImpl);
    mram_write(&desc, this->buffer + this->offset, sz);
    this->offsetInc(this, sz);
}

void xferItf_getUchar(xferItf *this, Uchar *val) { SCALAR_FETCH_MRAM_TO_WRAM(val, Uchar); }

void xferItf_getUllong(xferItf *this, Ullong *val) { SCALAR_FETCH_MRAM_TO_WRAM(val, Ullong); }
void xferItf_getInt(xferItf *this, Int *val) { SCALAR_FETCH_MRAM_TO_WRAM(val, Int); }

void xferItf_getVecUchar(xferItf *this, VecUchar *val, DescriptorImpl *len, enum DpuAlloctype atype)
{
    VECTOR_ALLOCATE_AND_FETCH_MRAM_TO_WRAM(val, len, Uchar, atype);
}
void xferItf_getVecChar(xferItf *this, char **val, DescriptorImpl *len, enum DpuAlloctype atype)
{
    VECTOR_ALLOCATE_AND_FETCH_MRAM_TO_WRAM(val, len, char, atype);
}

void xferItf_getVecCharNoAlloc(xferItf *this, DescriptorImpl *len)
{
    __dma_aligned DescriptorImpl type;
    this->getDesc(this, &type);
    this->getDesc(this, len);
    DescriptorImpl sz = XFER_MEM_ALIGN(*len * sizeof(char));
    this->offsetInc(this, sz);
}

void xferItf_getTensor2dUnsigned(xferItf *this, Tensor2d *t)
{
#define TYPE_BYTES_LOG2 2
#define TYPE_BYTES (1 << TYPE_BYTES_LOG2)

    __dma_aligned DescriptorImpl type;
    // type not checked yet
    this->getDesc(this, &type);
    this->getDesc(this, &(t->shapes[0]));
    this->getDesc(this, &(t->shapes[1]));
    t->aligned_shapes[0] = t->shapes[0];
    t->aligned_shapes[1] = XFER_MEM_ALIGN(t->shapes[1] << TYPE_BYTES_LOG2) >> TYPE_BYTES_LOG2;
    uint64_t aligned_size = t->aligned_shapes[0] * t->aligned_shapes[1];
    DescriptorImpl sz = XFER_MEM_ALIGN(aligned_size * TYPE_BYTES);
    t->mram_addr = buddy_alloc(t->shapes[0] * sizeof(uint8_t *));
    for (uint64_t i = 0; i < t->shapes[0]; i++) {
        t->mram_addr[i] = this->buffer + this->offset + (TYPE_BYTES * i * t->aligned_shapes[1]);
    }
    this->offsetInc(this, sz);

#undef TYPE_BYTES_LOG2
#undef TYPE_BYTES
}

/*
 * PUSH HAS A SPECIAL MEANING ON DPU side, it reserve the MRAM emplacement and solve the MRAM table of given tensor t
 * */
void xferItf_pushTensor2dUnsigned(xferItf *this, Tensor2d *t, uint64_t shapes[])
{
#define TYPE_BYTES_LOG2 2
#define TYPE_BYTES (1 << TYPE_BYTES_LOG2)

    // assign shapes first
    t->shapes[0] = shapes[0];
    t->shapes[1] = shapes[1];
    // type not checked yet
    this->pushDesc(this, CharType);
    this->pushDesc(this, (t->shapes[0]));
    this->pushDesc(this, (t->shapes[1]));
    t->aligned_shapes[1] = XFER_MEM_ALIGN(t->shapes[1] << TYPE_BYTES_LOG2) >> TYPE_BYTES_LOG2;
    uint64_t aligned_size = t->aligned_shapes[0] * t->aligned_shapes[1];
    DescriptorImpl sz = XFER_MEM_ALIGN(aligned_size * TYPE_BYTES);
    t->mram_addr = buddy_alloc(t->shapes[0] * sizeof(uint8_t *));
    for (uint64_t i = 0; i < t->shapes[0]; i++) {
        t->mram_addr[i] = this->buffer + this->offset + (TYPE_BYTES * i * t->aligned_shapes[1]);
    }
    this->offsetInc(this, sz);

#undef TYPE_BYTES_LOG2
#undef TYPE_BYTES
}

void xferItf_getTensor2dUchar(xferItf *this, Tensor2d *t)
{
#define TYPE_BYTES_LOG2 0
#define TYPE_BYTES (1 << TYPE_BYTES_LOG2)

    __dma_aligned DescriptorImpl type;
    // type not checked yet
    this->getDesc(this, &type);
    this->getDesc(this, &(t->shapes[0]));
    this->getDesc(this, &(t->shapes[1]));
    t->aligned_shapes[0] = t->shapes[0];
    t->aligned_shapes[1] = XFER_MEM_ALIGN(t->shapes[1] << TYPE_BYTES_LOG2) >> TYPE_BYTES_LOG2;
    uint64_t aligned_size = t->aligned_shapes[0] * t->aligned_shapes[1];
    DescriptorImpl sz = XFER_MEM_ALIGN(aligned_size * TYPE_BYTES);
    t->mram_addr = buddy_alloc(t->shapes[0] * sizeof(uint8_t *));
    for (uint64_t i = 0; i < t->shapes[0]; i++) {
        t->mram_addr[i] = this->buffer + this->offset + (TYPE_BYTES * i * t->aligned_shapes[1]);
    }
    this->offsetInc(this, sz);

#undef TYPE_BYTES_LOG2
#undef TYPE_BYTES
}

/*
 * PUSH HAS A SPECIAL MEANING ON DPU side, it reserve the MRAM emplacement and solve the MRAM table of given tensor t
 * */
void xferItf_pushTensor2dUchar(xferItf *this, Tensor2d *t, uint64_t shapes[])
{
#define TYPE_BYTES_LOG2 0
#define TYPE_BYTES (1 << TYPE_BYTES_LOG2)

    // assign shapes first
    t->shapes[0] = shapes[0];
    t->shapes[1] = shapes[1];
    // type not checked yet
    this->pushDesc(this, CharType);
    this->pushDesc(this, (t->shapes[0]));
    this->pushDesc(this, (t->shapes[1]));
    t->aligned_shapes[0] = t->shapes[0];
    t->aligned_shapes[1] = XFER_MEM_ALIGN(t->shapes[1] << TYPE_BYTES_LOG2) >> TYPE_BYTES_LOG2;
    uint64_t aligned_size = t->aligned_shapes[0] * t->aligned_shapes[1];
    DescriptorImpl sz = XFER_MEM_ALIGN(aligned_size * TYPE_BYTES);
    t->mram_addr = buddy_alloc(t->shapes[0] * sizeof(uint8_t *));
    for (uint64_t i = 0; i < t->shapes[0]; i++) {
        t->mram_addr[i] = this->buffer + this->offset + (TYPE_BYTES * i * t->aligned_shapes[1]);
    }
    this->offsetInc(this, sz);

#undef TYPE_BYTES_LOG2
#undef TYPE_BYTES
}

void xferItf_pushTensor2dUINT64(xferItf *this, Tensor2d *t, uint64_t shapes[])
{
#define TYPE_BYTES_LOG2 3
#define TYPE_BYTES (1 << TYPE_BYTES_LOG2)

    // assign shapes first
    t->shapes[0] = shapes[0];
    t->shapes[1] = shapes[1];
    // type not checked yet
    this->pushDesc(this, CharType);
    this->pushDesc(this, (t->shapes[0]));
    this->pushDesc(this, (t->shapes[1]));
    t->aligned_shapes[0] = t->shapes[0];
    t->aligned_shapes[1] = XFER_MEM_ALIGN(t->shapes[1] << TYPE_BYTES_LOG2) >> TYPE_BYTES_LOG2;
    uint64_t aligned_size = t->aligned_shapes[0] * t->aligned_shapes[1];
    DescriptorImpl sz = XFER_MEM_ALIGN(aligned_size * TYPE_BYTES);
    t->mram_addr = buddy_alloc(t->shapes[0] * sizeof(uint8_t *));
    for (uint64_t i = 0; i < t->shapes[0]; i++) {
        t->mram_addr[i] = this->buffer + this->offset + (TYPE_BYTES * i * t->aligned_shapes[1]);
    }
    this->offsetInc(this, sz);

#undef TYPE_BYTES_LOG2
#undef TYPE_BYTES
}

void xferItf_getVecInt(xferItf *this, VecInt *val, DescriptorImpl *len, enum DpuAlloctype atype)
{
    VECTOR_ALLOCATE_AND_FETCH_MRAM_TO_WRAM(val, len, Int, atype);
}
void xferItf_getVecUllong(xferItf *this, VecUllong *val, DescriptorImpl *len, enum DpuAlloctype atype)
{
    VECTOR_ALLOCATE_AND_FETCH_MRAM_TO_WRAM(val, len, Ullong, atype);
}
void xferItf_pushVecUchar(xferItf *this, VecUchar val, DescriptorImpl len)
{
    VECTOR_FETCH_WRAM_TO_MRAM(val, len, Uchar, UcharType);
}

int64_t DpuAllign(uint64_t val) { return XFER_MEM_ALIGN(val); }

#undef VECTOR_ALLOCATE_AND_FETCH_MRAM_TO_WRAM
#undef VECTOR_FETCH_WRAM_TO_MRAM
#undef SCALAR_FETCH_MRAM_TO_WRAM

#endif

#endif /* -_XFERITF_- */
