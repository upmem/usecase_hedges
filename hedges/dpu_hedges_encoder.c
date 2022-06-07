/**
 * @file dpu_hedges_encoder.c
 * @author Dimitri Gerin (dgerin@upmem.com)
 * @brief DPU implementation of HEDGES inner encoder
 *        # original project #
 *        Based on original HEDGES project
 *        HEDGES Error-Correcting Code for DNA Storage Corrects Indels and
 *        Allows Sequence Constraints.
 *        William H. Press, John A. Hawkins, Stephen Knox Jones Jr,
 *        Jeffrey M. Schaub, Ilya J. Finkelstein
 *        submitted to Proceedings of the National Academy of Sciences.
 * @copyright 2022 UPMEM
 */

#include <barrier.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>
#include <stdint.h>
#include <stdio.h>
#include <alloc.h>
#include <assert.h>
#include <perfcounter.h>
#include <seqread.h>
#include <common.h>
#include <../hedges/ran.h>
#include <xferItf.h>

uint8_t __mram_noinit xferitf_buffer[1 << HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE];
__host uint64_t loaded = 0;
__host uint64_t xferOffset;

/* HEDGES Encoder constant parameters */
__host Int HLIMIT;
__host double _reward = -0.13;
__host Int NSTAK;

__host uint32_t NPACKET_PER_DPU;

/**
 * MAXSEQ_DPU is used to set buffer size of both unpacked
 * bits buffer, and uncoded bit buffer. Theses two buffer
 * are declared in TASKLET stack, so we must limit its size
 * NOTE that the effective size of both buffer are computed as
 * uint64_t codetext_size = vbitlen(8 / 8 bit per packed bit / packed_message_size) + RPRIMER;
 * but we must use constants for C array SIZE, so we fix the limit
 * at the smallest possible value right now.
 * WARNING : 300 works only for coderate = 0.5
 */

#define MAXSEQ_DPU ENCODER_STACK_MAXSEQ_DPU_NBYTES

__dma_aligned Int NPREV;
__dma_aligned Int HSALT;
__dma_aligned Int LPRIMER;
__dma_aligned Int RPRIMER;
__dma_aligned Ullong prevmask;
__dma_aligned Ullong seqnomask;
__dma_aligned Ullong saltmask;
__dma_aligned Int NSP;
__dma_aligned Int DNAWINDOW;
__dma_aligned Int MAXGC;
__dma_aligned Int MINGC;
__dma_aligned Int MAXRUN;
__dma_aligned GF4reg dnawinmask;
__dma_aligned GF4reg dnaoldmask;
__dma_aligned GF4reg acgtacgt;

__dma_aligned uint64_t RIGHTPRIMER_LEN;
__dma_aligned uint64_t PRIMERSALT_LEN;
__dma_aligned uint64_t PATTARR_LEN;

__dma_aligned GF4word rightprimer;
__dma_aligned VecUllong primersalt;
__dma_aligned VecInt pattarr;

/* input tensor */
Tensor2d I;
/* output tensor */
Tensor2d O;

/**
 * @brief returns the number of allowed ACGTs and puts them in dnac_ok
 *
 * @param dnac_ok
 * @param prev
 * @return returns
 */
Int dnacallowed(Uchar *dnac_ok, GF4reg prev)
{
    /* returns the number of allowed ACGTs and puts them in dnac_ok */
    if (DNAWINDOW <= 0) {
        dnac_ok[0] = 0;
        dnac_ok[1] = 1;
        dnac_ok[2] = 2;
        dnac_ok[3] = 3;
        return 4;
    }
    Int ans, gccount, last = prev & 3, nrun = 1;
    bool isrun = false;
    Ullong reg;
    /* get GCcount */
    reg = prev & dnaoldmask;
    /* makes ones for GC, zeros for AT */
    reg = (reg ^ (reg >> 1)) & 0x5555555555555555ull;
    /* popcount inline */
    reg -= ((reg >> 1) & 0x5555555555555555ull);
    reg = (reg & 0x3333333333333333ull) + (reg >> 2 & 0x3333333333333333ull);
    /* the popcount */
    gccount = ((reg + (reg >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
    /* is there a run and, if so, of what */
    reg = (prev >> 2);
    while ((reg & 3) == last) {
        ++nrun;
        if (nrun >= MAXRUN) {
            isrun = true;
            break;
        }
        reg >>= 2;
    }
    /* the horrible logic tree */
    if (gccount >= MAXGC) {
        ans = 2;
        /* A is ok */
        dnac_ok[0] = 0;
        /* T is ok */
        dnac_ok[1] = 3;
        if (isrun) {
            if (last == 0) {
                ans = 1;
                /* only T ok */
                dnac_ok[0] = 3;
            } else if (last == 3) {
                ans = 1;
                /* only A ok */
                dnac_ok[0] = 0;
            }
        }
    } else if (gccount <= MINGC) {
        ans = 2;
        /* C is ok*/
        dnac_ok[0] = 1;
        dnac_ok[1] = 2;
        /* G is ok*/
        if (isrun) {
            if (last == 1) {
                ans = 1;
                /* only G ok*/
                dnac_ok[0] = 2;
            } else if (last == 2) {
                ans = 1;
                /* only C ok*/
                dnac_ok[0] = 1;
            }
        }
    } else { /* no GC constraints*/
        ans = 4;
        /* A is ok */
        dnac_ok[0] = 0;
        /* C is ok */
        dnac_ok[1] = 1;
        /* G is ok */
        dnac_ok[2] = 2;
        /* T is ok */
        dnac_ok[3] = 3;
        if (isrun) {
            ans = 3;
            for (int i = last; i < 3; i++)
                dnac_ok[i] = dnac_ok[i + 1];
        }
    }
    return ans;
}

/**
 * @brief HEDGEES INNER encoder/decoder hash function
 *
 * @param out

 * @param bits
 * @param seq
 * @param salt
 * @param mod
 */
void hash_dna(Int *out, Ullong bits, Int seq, Ullong salt, Int mod)
{
    Ullong ret;

    ranhash_int64(&ret, ((((((Ullong)(seq)&seqnomask) << NPREV) | bits) << HSALT) | salt));

    Ullong div;
    Ullong ret_;

    /* op: (% mod) with mod from 0 to 4 */
    if (mod == 4) {
        div = ret >> 2;
        ret_ = ret - (div << 2);
    } else if (mod == 2) {
        div = ret >> 1;
        ret_ = ret - (div << 1);
    } else if (mod == 1) {
        ret_ = ret;
    } else {
        ret_ = ret % mod;
    }

    *out = (Int)(ret_);
}

/**
 * @brief  Computes the message len from vbits mapping patern.
 * @param nmb
 * @return Int
 */
Int vbitlen(Int nmb)
{ /* how long is message in vbits?  (patarr must already be set) */
    Int ksize, nn = 0;
    for (ksize = 0;; ksize++) { /* how many Mbits do we need? */
        if (nn >= nmb)
            break;
        assert(ksize < ENCODER_STACK_MAXSEQ_DPU_NBYTES && "vbitlen: MAXSEQ too small");
        nn += pattarr[ksize];
    }
    return ksize;
}

/**
 * @brief transform efficiently encoded bit sequence into sequence of bit stored as byte
 *
 * @param packed_bytes
 * @param ans output bit sequence
 * @param message input message
 * @param n message size
 */
void unpackvbits(uint64_t *unpacked_len, Uchar *ans, Uchar *message, Int n)
{
    Int i, j, nmb = 8 * n, k, k1;
    uint64_t ksize;
    Uchar bit;
    // TODO this function could be called only once at programm start
    ksize = vbitlen(nmb);

    for (k = 0; k < ksize; k++)
        ans[k] = 0;

    i = j = 0;
    for (k = 0; k < ksize; k++) {
        for (k1 = 0; k1 < pattarr[k]; k1++) {
            bit = (i < n ? (message[i] >> (7 - j++)) & 1 : 0);
            if (j == 8) {
                j = 0;
                ++i;
            }
            ans[k] = (ans[k] << 1) | bit;
        }
    }
    *unpacked_len = ksize;
}

/**
 * @brief encode one input sequence (from binary space to ACGT space)
 *
 * @param packet_index_ absolute index of the input sequence to encode into the input packet
 * @return
 */
void encode(uint32_t packet_index_)
{
    Int regout;

    Uchar vbits[MAXSEQ_DPU];
    __dma_aligned Uchar message[MAXSEQ_DPU];
    uint64_t message_len = I.shapes[1];

    /*
     * get one line of Input 2d tensor that is
     * processed by one tasklet
     * */

    uint64_t input_line_size = sizeof(char) * I.aligned_shapes[1];
    mram_read(I.mram_addr[packet_index_], message, input_line_size);
    /*
     * unpacked sequence len, computed by unpackvbits
     * (variable param to support different code rate)
     * */
    uint64_t unpacked_message_size;

    unpackvbits(&unpacked_message_size, vbits, message, message_len);
    Int k = 0, nbits, mod;

    uint64_t output_line_size = sizeof(Uchar) * O.aligned_shapes[1];
    __dma_aligned Uchar codetext[MAXSEQ_DPU];

    Mbit messagebit;
    Ullong prevbits = 0, salt = 0, newsalt = 0;
    /* initialize with no runs and balanced cg */
    GF4reg prevcode = acgtacgt;

    Uchar dnac_ok[4];
    for (k = 0; k < unpacked_message_size; k++) { /* on decoding, k is called seq */
        messagebit = vbits[k];
        nbits = pattarr[k];
        if (k < LPRIMER) {
            salt = primersalt[k];
        } else if (k < NSP) {
            salt = 0;
            newsalt = ((newsalt << 1) & saltmask) ^ messagebit;
        } else if (k == NSP) {
            /* time to update the salt */
            salt = newsalt;
        }
        mod = (k < LPRIMER ? 4 : dnacallowed(dnac_ok, prevcode));
        hash_dna(&regout, prevbits, k, salt, mod);
        regout = (regout + (Uchar)(messagebit)) % mod;

        codetext[k] = (k < LPRIMER ? regout : dnac_ok[regout]);
        prevcode = ((prevcode << 2) | codetext[k]) & dnawinmask;
        /* variable number */
        prevbits = ((prevbits << nbits) & prevmask) | messagebit;
    }
    for (k = 0; k < RPRIMER; k++) {
        codetext[k + unpacked_message_size] = rightprimer[k];
    }
    mram_write(codetext, O.mram_addr[packet_index_], output_line_size);
}

BARRIER_INIT(barrier, NR_TASKLETS);
uint32_t packet_index = 0;
MUTEX_INIT(packet_index_lock);

/**
 * @brief concurently provide next packet index
 *
 * @param index
 * @return
 */
void get_next_packet_index(uint32_t *index)
{
    mutex_lock(packet_index_lock);
    *index = packet_index;
    packet_index++;
    mutex_unlock(packet_index_lock);
}

/**
 * @brief reset packet index
 *
 * @param limit
 * @param nmessbits
 * @param input_codetext
 * @return
 */
void reset_packet_index()
{
    mutex_lock(packet_index_lock);
    packet_index = 0;
    mutex_unlock(packet_index_lock);
}
/* declare HOST/DPU xfer interface */
xferItf xitf = XFERITF_INIT();

/**
 * @brief DPU HEDGES encoder entry point
 * @param nmessbits number of bit in each sequence
 * @param global I -> input packet of sequence
 *
 * @return
 */
int main()
{
    if (me() == 0) {
        xitf.init(&xitf, xferitf_buffer);

        if (!loaded) {
            buddy_init(BUDDY_ALLOCATOR_SIZE_BYTE);

            xitf.getInt(&xitf, &NPREV);
            xitf.getInt(&xitf, &HSALT);
            xitf.getInt(&xitf, &LPRIMER);
            xitf.getInt(&xitf, &RPRIMER);

            xitf.getUllong(&xitf, &prevmask);
            xitf.getUllong(&xitf, &seqnomask);
            xitf.getUllong(&xitf, &saltmask);

            xitf.getInt(&xitf, &NSP);
            xitf.getInt(&xitf, &DNAWINDOW);
            xitf.getInt(&xitf, &MAXGC);
            xitf.getInt(&xitf, &MINGC);
            xitf.getInt(&xitf, &MAXRUN);

            xitf.getUllong(&xitf, &dnawinmask);
            xitf.getUllong(&xitf, &dnaoldmask);
            xitf.getUllong(&xitf, &acgtacgt);
            xitf.getVecUchar(&xitf, &rightprimer, &RIGHTPRIMER_LEN, mem);
            xitf.getVecUllong(&xitf, &primersalt, &PRIMERSALT_LEN, mem);
            xitf.getVecInt(&xitf, &pattarr, &PATTARR_LEN, mem);
            xitf.save(&xitf, &xferOffset);

            xitf.restore(&xitf, &xferOffset);
            reset_packet_index();

            xitf.getTensor2dUchar(&xitf, &I);

            /* packed bits are stored as rows */
            uint64_t packed_message_size = I.shapes[1];
            /* number of encoded output character (in output alphabet space) */
            uint64_t codetext_size = vbitlen(8 /* 8 bit per packet */ * packed_message_size) + RPRIMER;
            uint64_t output_shapes[2] = { I.shapes[0], codetext_size };

            /* because push() in DPU side, push means : allocate MRAM space + solve mram table */
            xitf.pushTensor2dUchar(&xitf, &O, output_shapes);
            loaded = 1;
        }
    }
    barrier_wait(&barrier);

    printf("tasklet %u: stack = %u\n", me(), check_stack());
    uint32_t tasklet_current_packet_index;
    uint32_t total_packet = O.shapes[0];

    /* call encode for each packet */
    while (1) {
        get_next_packet_index(&tasklet_current_packet_index);
        if (tasklet_current_packet_index >= total_packet)
            break;
        encode(tasklet_current_packet_index);
    }

    barrier_wait(&barrier);

    if (me() == 0) {
        free2dTensor(&I);
        free2dTensor(&O);
        reset_packet_index();

        printf("END\n");
    }
    return 0;
}
