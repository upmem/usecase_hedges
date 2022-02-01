#include <barrier.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>
#include <stdint.h>
#include <stdio.h>
#include <alloc.h>
#include <assert.h>
#include <seqread.h>
#include <common.h>
#include <../hedges/ran.h>
#include <xferItf.h>
#include <profiling.h>
#include "heap_1_1.h"

#include <perfcounter.h>

#if defined(MESURE_PERF)
__host perfcounter_t nb_cycles_hypcompute[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_decode[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_io[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_hashfunc[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_penality[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_hypcompute_first_section[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_hypload[NR_TASKLETS] = {0};

//__host perfcounter_t hypcompute_time = 0;
//__host perfcounter_t hypcompute_prct;
#endif

#ifdef MESURE_BW
__host uint64_t nb_bytes_loaded[NR_TASKLETS] = {0};
__host uint64_t nb_bytes_written[NR_TASKLETS] = {0};
#endif

typedef struct Hypothesis
{
  /** next char in message **/
  Int offset;
  /** my position in the decoded message (0,1,...) **/
  Int seq;
  /** last decoded up to now **/
  Mbit messagebit;
  /** corresponding salt values **/
  Ullong prevbits, salt, newsalt;
  /** precedent code **/
  GF4reg prevcode;
  /** cumulated score **/
  heap_score_type score;
  /** index of predecessor in hypostackp **/
  Int predi;
  uint64_t tmp;
} hypothesis;

__mram_noinit uint8_t hstack[NR_TASKLETS][ XFER_MEM_ALIGN(sizeof(hypothesis)  *HEAP_MAX_ITEM  )];
__dma_aligned uint64_t heap_pos[NR_TASKLETS];



__dma_aligned uint64_t heap_pos_final[NR_TASKLETS];
__mram_noinit uint8_t xferitf_buffer[1 << HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE];
__host __dma_aligned uint64_t loaded = 0;
__host __dma_aligned uint64_t xferOffset;
__host uint64_t nr_tasklets = NR_TASKLETS;

//__mram_noinit uint8_t scores[NR_TASKLETS * XFER_MEM_ALIGN(sizeof(heap_score_type) * HEAP_MAX_ITEM)];
//__mram_noinit uint8_t ptrs[NR_TASKLETS * XFER_MEM_ALIGN(sizeof(heap_ptr_type) * HEAP_MAX_ITEM)];
__mram_noinit uint8_t heap_items[NR_TASKLETS * XFER_MEM_ALIGN(sizeof(item) * HEAP_MAX_ITEM)];

__dma_aligned uint64_t nhypo[NR_TASKLETS];
__dma_aligned heap heap_;

// MAXSEQ_DPU is used to set buffer size of both unpacked
// bits buffer, and uncoded bit buffer. Theses two buffer
// are declared in TASKLET stack, so we must limit its size
// NOTE that the effective size of both buffer are computed as
// uint64_t codetext_size = vbitlen(8 /*8 bit per packed bit*/ * packed_message_size) + RPRIMER;
// but we must use constants for C array SIZE, so we fix the limit
// at the smallest possible value right now.
// WARNING : 300 works only for coderate = 0.5
#define MAXSEQ_DPU ENCODER_STACK_MAXSEQ_DPU_NBYTES

__host uint64_t total_cycles[NR_TASKLETS];
__host uint64_t perf_type;
__host uint64_t perf_bw;
// __dma_aligned Uchar message[NR_TASKLETS][MAXSEQ_DPU];

/**
 * HOST/DPU shared global variables
 * **/
/**
 * HEDGES Decoder constant parameters
 * **/
__host Int HLIMIT;
__host Int NSTAK;

__host uint64_t nmessbit;

// constants seted at programm start
Int codetextlen_g;

// these are the rewards and penalties applied at each position
__host double _reward = -0.13;
double _substitution = 1.;
double _deletion = 1.;
double _insertion = 1.;
double _dither = 0.;

#define DOUBLE_TO_FP(i) (heap_score_type)((double)(i) * ((double)(1 << DECODER_QUANT_FRAC_BITS)))

heap_score_type reward_;
heap_score_type substitution_;
heap_score_type deletion_;
heap_score_type insertion_;
heap_score_type dither_;

void quantif_test()
{
  double a = 0.0;
  double b;
  b = a;
  for (uint64_t i = 0; i < 10; i++)
  {
    b = b + _reward;
    printf("b : %lf\n", b);
  }
}

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

Tensor2d I;
Tensor2d O;
Tensor2d O_scores, O_ptr;

Int dnacallowed(Uchar *dnac_ok, GF4reg prev)
{
  // returns the number of allowed ACGTs and puts them in dnac_ok
  if (DNAWINDOW <= 0)
  {
    dnac_ok[0] = 0;
    dnac_ok[1] = 1;
    dnac_ok[2] = 2;
    dnac_ok[3] = 3;
    return 4;
  }
  Int ans, gccount, last = prev & 3, nrun = 1;
  bool isrun = false;
  Ullong reg;
  // get GCcount
  reg = prev & dnaoldmask;
  reg = (reg ^ (reg >> 1)) & 0x5555555555555555ull; // makes ones for GC, zeros for AT
  // popcount inline:
  reg -= ((reg >> 1) & 0x5555555555555555ull);
  reg = (reg & 0x3333333333333333ull) + (reg >> 2 & 0x3333333333333333ull);
  gccount = ((reg + (reg >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56; // the popcount
  // is there a run and, if so, of what
  reg = (prev >> 2);
  while ((reg & 3) == last)
  {
    ++nrun;
    if (nrun >= MAXRUN)
    {
      isrun = true;
      break;
    }
    reg >>= 2;
  }
  // the horrible logic tree:
  if (gccount >= MAXGC)
  {
    ans = 2;
    dnac_ok[0] = 0; // A is ok
    dnac_ok[1] = 3; // T is ok
    if (isrun)
    {
      if (last == 0)
      {
        ans = 1;
        dnac_ok[0] = 3; // only T ok
      }
      else if (last == 3)
      {
        ans = 1;
        dnac_ok[0] = 0; // only A ok
      }
    }
  }
  else if (gccount <= MINGC)
  {
    ans = 2;
    dnac_ok[0] = 1; // C is ok
    dnac_ok[1] = 2; // G is ok
    if (isrun)
    {
      if (last == 1)
      {
        ans = 1;
        dnac_ok[0] = 2; // only G ok
      }
      else if (last == 2)
      {
        ans = 1;
        dnac_ok[0] = 1; // only C ok
      }
    }
  }
  else
  { // no GC constraints
    ans = 4;
    dnac_ok[0] = 0; // A is ok
    dnac_ok[1] = 1; // C is ok
    dnac_ok[2] = 2; // G is ok
    dnac_ok[3] = 3; // T is ok
    if (isrun)
    {
      ans = 3;
      for (int i = last; i < 3; i++)
        dnac_ok[i] = dnac_ok[i + 1];
    }
  }
  return ans;
}

Int digest(Ullong bits, Int seq, Ullong salt, Int mod)
{
#if defined(MESURE_PERF)
  perfcounter_t start_time;
  if (perf_type)
    start_time = perfcounter_get();
#endif

  Int out = (Int)(ranhash_int64((((((Ullong)(seq)&seqnomask) << NPREV) | bits) << HSALT) | salt) % mod);

#if defined(MESURE_PERF)
  if (perf_type)
    nb_cycles_hashfunc[me()] += perfcounter_get() - start_time;
#endif

  return out;
}

/**
 * Computes the message len from vbits mapping patern.
 * TODO : clean patarr global var MAXSEQ_DPU
 * **/
Int vbitlen(Int nmb)
{
  Int ksize, nn = 0;
  for (ksize = 0;; ksize++)
  {
    if (nn >= nmb)
      break;
    assert(ksize < MAXSEQ_DPU && "vbitlen: MAXSEQ too small");
    nn += pattarr[ksize];
  }
  return ksize;
}

void packvbits(Uchar *packed_bytes, Uchar *vbits, uint64_t nmessbits, uint64_t vbits_size)
{
  Int i, j, k, k1;
  Uchar bit;
  uint64_t nn;
  assert(vbits_size <= MAX_DECODED_VBITS);
  for (k = 0; k < vbits_size; k++)
    nn += pattarr[k];      // number of bits
  nn = MIN(nn, nmessbits); // no more than the specified number of bits
  nn = nmessbits;
  uint64_t nbytes = (nn + 7) / 8; // number of bytes
  // printf(" nbytes %lu, vbits size %lu \n", nbytes, vbits_size);

  i = j = 0;
  for (k = 0; k < vbits_size; k++)
  {
    for (k1 = pattarr[k] - 1; k1 >= 0; k1--)
    {
      bit = (vbits[k] >> k1) & 1;
      packed_bytes[i] = packed_bytes[i] | (bit << (7 - j++));
      if (j == 8)
      {
        j = 0;
        if (++i == nbytes)
          break;
      }
    }
    if (i == nbytes)
      break;
  }
}

xferItf xitf = XFERITF_INIT();

#define MAX_CODETEXT_LEN 300

/**
 * cache for decoded input sequence
 **/
//__dma_aligned Uchar input_codetext[NR_TASKLETS][MAX_CODETEXT_LEN];

void hypothesis_init_root(hypothesis *h)
{
  h->predi = -1;
  h->offset = -1;
  h->seq = -1;
  h->messagebit = 0;
  h->prevbits = 0;
  h->score = 0.;
  h->salt = 0;
  h->newsalt = 0;
  h->prevcode = acgtacgt;
}


void hypothesis_load(__dma_aligned hypothesis *cache,  uint32_t pos)
{
  //printf("r %u \n", pos);


  
  mram_read(&(hstack[me()][pos * sizeof(hypothesis) ]), cache, XFER_MEM_ALIGN(sizeof(hypothesis)));
#ifdef MESURE_BW
  if (perf_bw)
    nb_bytes_loaded[me()] += sizeof(hypothesis);
#endif
}

void hypothesis_store(__dma_aligned hypothesis *cache, uint32_t pos)
{

  //printf("w %u \n", pos);
  mram_write(cache, &(hstack[me()][pos * sizeof(hypothesis)]), XFER_MEM_ALIGN(sizeof(hypothesis)));
#ifdef MESURE_BW
  if (perf_bw)
    nb_bytes_written[me()] += sizeof(hypothesis);
#endif
}

// void init_heap_and_stack()
// {
//  __dma_aligned hypothesis cached_h;
// 
//   heap_init();
// 
// hypothesis_init_root(&cached_h);
//  hypothesis_store(&cached_h, me(), 0);
// 
//   /** set nhypo to 1 element **/
//   nhypo[me()] = 1;
// 
//   heap_pos[me()] = 0;
// heap_pos_final[me()] = 0;
// 
//   __dma_aligned heap_score_type score = 10000000000;
//   __dma_aligned heap_ptr_type ptr = 0;
//   heap_push(&score, &ptr, perf_type);
//   return;
// }

void reset_heap_and_stack()
{
  __dma_aligned heap_ptr_type zero = 0;
  __dma_aligned hypothesis cached_h;
  __dma_aligned item cur_heap_item;

  hypothesis_init_root(&cached_h);
  hypothesis_store(&cached_h,  0);

  /** set nhypo to 1 element **/
  nhypo[me()] = 1;
  heap_.heap_pos[me()] = 0;



  //for (uint64_t i = 0; i < HEAP_MAX_ITEM; i++) {
  //  heap_write_score(&heap_, i, &score_max_val);
  //  heap_write_ptr(&heap_, i, &zero);
  //}


  cur_heap_item.score = 1000000000; // 1410065408
  cur_heap_item.ptr = 0;
  heap_push(&heap_, &cur_heap_item, perf_type);


  return;
}

Int init_from_predecessor(hypothesis *h, Int pred, Mbit mbit, Int skew, __dma_aligned Uchar *input_codetext)
{
#if defined(MESURE_PERF)
  perfcounter_t start_time;
  if (perf_type)
    start_time = perfcounter_get();
#endif

  bool discrep;
  Int regout, mod;
  heap_score_type mypenalty;
  Ullong mysalt;

#if defined(MESURE_PERF)
  perfcounter_t start_time__;
  if (perf_type)
    start_time__ = perfcounter_get();
#endif

  /** load precedent hypothesis**/
  __dma_aligned hypothesis hp;
  hypothesis_load(&hp,  pred);

#if defined(MESURE_PERF)
  if (perf_type)
    nb_cycles_hypload[me()] += perfcounter_get() - start_time__;
#endif

#if defined(MESURE_PERF)
  perfcounter_t start_time___;
  if (perf_type)
    start_time___ = perfcounter_get();
#endif

  Uchar dnac_ok[4];
  h->predi = pred;
  h->messagebit = mbit; // variable number
  h->seq = hp.seq + 1;
  if (h->seq > MAXSEQ_DPU)
    assert(false && "init_from_predecessor: MAXSEQ too small");
  Int nbits = pattarr[h->seq];
  h->prevbits = hp.prevbits;
  h->salt = hp.salt;
  if (h->seq < LPRIMER)
  {
    mysalt = primersalt[h->seq];
  }
  else if (h->seq < NSP)
  {
    mysalt = h->salt;
    h->newsalt = ((hp.newsalt << 1) & saltmask) ^ h->messagebit; // variable bits overlap, but that's ok with XOR
  }
  else if (h->seq == NSP)
  {
    /** update salt value **/
    mysalt = h->salt = hp.newsalt;
  }
  else
    mysalt = h->salt;
  h->offset = hp.offset + 1 + skew;
  if (h->offset >= codetextlen_g)
    return 0;

  /** calculate predicted message under this hypothesis **/
  h->prevcode = hp.prevcode;
  mod = (h->seq < LPRIMER ? 4 : dnacallowed(dnac_ok, h->prevcode));

#if defined(MESURE_PERF)
  if (perf_type)
    nb_cycles_hypcompute_first_section[me()] += perfcounter_get() - start_time___;
#endif

  regout = digest(h->prevbits, h->seq, mysalt, mod);
  regout = (regout + (Uchar)(h->messagebit)) % mod;
  regout = (h->seq < LPRIMER ? regout : dnac_ok[regout]);
  h->prevbits = ((hp.prevbits << nbits) & prevmask) | h->messagebit; // variable number
  h->prevcode = ((h->prevcode << 2) | regout) & dnawinmask;

#if defined(MESURE_PERF)
  perfcounter_t start_time_;
  if (perf_type)
    start_time_ = perfcounter_get();
#endif

  /**deletion**/
  if (skew < 0)
  {
    mypenalty = deletion_;
  }
  else
  {
    /** the only place where a check is possible!  **/
    discrep = (regout == input_codetext[h->offset]);
    /** substitution of error free**/
    if (skew == 0)
      mypenalty = (discrep ? reward_ : substitution_);
    else
    {
      /** inserion **/
      mypenalty = insertion_ + (discrep ? reward_ : substitution_);
    }
  }

#if defined(MESURE_PERF)
  if (perf_type)
    nb_cycles_penality[me()] += perfcounter_get() - start_time_;
#endif

  /** NOTE : dont't understant why this code does.
   * It adds a randomly signed zero -0./+O.
   * TOTO : understand the impact of this step
   * in the scoring.
   *
   * if (dither > 0.)
   *  mypenalty += dither * (2. * ran_double() - 1.);
   **/
  h->score = hp.score + mypenalty;

#if defined(MESURE_PERF)
  if (perf_type)
    nb_cycles_hypcompute[me()] += perfcounter_get() - start_time;
#endif

  return 1;
}

void shoveltheheap(Int limit, Int nmessbits, __dma_aligned Uchar *input_codetext)
{
  // given the heap, keep processing it until offset limit, hypothesis limit, or an error is reached
  Int seq, nguess, qqmax = -1, ofmax = -1, seqmax = vbitlen(nmessbits);
  Uchar mbit;
  __dma_aligned item cur_heap_item;
  __dma_aligned hypothesis hp;
  __dma_aligned hypothesis new_hp;
  __dma_aligned item new_heap_item;
  uint8_t errcode = 0;
  cur_heap_item.score = 0;
 cur_heap_item.ptr = 0;

  while (true)
  {
    /**
     * pop heap -> gives best curscore and its associated
     * hypothesis position (hypothesis stack)
     * **/
    bool nemptyheap = heap_pop(&heap_, &cur_heap_item, perf_type);

    assert(nemptyheap && "pop empty heap");

    /**load hypothesis from stack**/
    hypothesis_load(&hp,  cur_heap_item.ptr);

    seq = hp.seq;
    assert(seq < MAXSEQ_DPU && "shoveltheheap: MAXSEQ too small");

    nguess = 1 << pattarr[seq + 1]; // i.e., 1, 2, or 4
    if (hp.offset > ofmax)
    { // keep track of farthest gotten to
      ofmax = hp.offset;
      qqmax = cur_heap_item.ptr;
    }
    if (cur_heap_item.score > 1000000000)
      break; // heap is empty
    if (hp.offset >= limit - 1)
      break; // errcode 0 (nominal success)
    if (nmessbits > 0 && seq >= seqmax - 1)
      break; // ditto when no. of message bits specified

    // nhypo below gives a new position in the HEAP
    if (nhypo[me()] > HLIMIT)
    {
      errcode = 2;
      heap_pos_final[me()] = qqmax;
      return;
    }
    /**
     *   NOTE : heap/stack resising, this feature doesn t exists on DPU
     *   if (nhypo + 12 >= nnstak)
     *  {
     *    nnstak *= 2;
     *    hypostack.resize(nnstak, true);
     *    if (hypostack.size() != nnstak)
     *      NRpyException("resize of hypostack failed");
     *  }
     **/

    /** sequential computing of each kind of Hypothesis, for each bits **/
    /** error free hypothesis **/
    for (mbit = 0; mbit < nguess; mbit++)
    {
      if (init_from_predecessor(&new_hp, cur_heap_item.ptr , mbit, 0, input_codetext))
      {
        new_heap_item.score = new_hp.score;
        new_heap_item.ptr = nhypo[me()];

        heap_push(&heap_, &new_heap_item, perf_type);
        hypothesis_store(&new_hp , nhypo[me()]);
        nhypo[me()]++;
      }
    }
    /** deletion error hypothesis **/
    for (mbit = 0; mbit < nguess; mbit++)
    {
      if (init_from_predecessor(&new_hp, cur_heap_item.ptr , mbit, -1, input_codetext))
      {
        new_heap_item.score = new_hp.score;
        new_heap_item.ptr = nhypo[me()];

        heap_push(&heap_, &new_heap_item, perf_type);
        hypothesis_store(&new_hp, nhypo[me()]);
        nhypo[me()]++;
      }
    }
    /** insertion error hypothesis **/
    for (mbit = 0; mbit < nguess; mbit++)
    {
      if (init_from_predecessor(&new_hp, cur_heap_item.ptr , mbit, 1, input_codetext))
      {
        new_heap_item.score = new_hp.score;
        new_heap_item.ptr = nhypo[me()];

        heap_push(&heap_, &new_heap_item, perf_type);
        hypothesis_store(&new_hp, nhypo[me()]);
        nhypo[me()]++;
      }
    }
  }
  heap_pos_final[me()] = cur_heap_item.ptr;
}

void print_heap_final_pos()
{
  printf("nhyp[%u] %lu \n", (unsigned)me(), heap_pos_final[me()]);
}

/**
 * extract vbits from heap
 * **/
void traceback(__dma_aligned Uchar *vbits, uint64_t *traceback_size)
{

  Int k, kk = 0;
  Int nfinal = (Int)(heap_pos_final[me()]);
  Int q = nfinal;
  // printf("final in traceback [%d] %d \n",  me(), q);
  __dma_aligned hypothesis h;
  /**
   * Computes the length of the decoded bits sequence
   * starting from the heap final position.
   * **/
  Int q_;
  q_ = q;
  do
  {
    hypothesis_load(&h,  q);

    //printf("final in traceback [%d] %d  h.perdi %d score %d \n", me(), q,  h.predi , h.score  );
    if (!(h.predi > 0))
      break;
    kk++;
    //if (kk>4)
     // return;
    q_ = q;
    q = h.predi;
    // if (q_ == q)
    //   return;
  } while (h.predi > 0);

  *traceback_size = kk + 1;
  q = nfinal;
  k = kk;

  /**
   * Copy each decoded bit of decoded sequence
   * from hypothesis stack to output array (vbits)
   **/
  do
  {
    hypothesis_load(&h,  q);
    if (!(h.predi > 0))
      break;
    vbits[k] = h.messagebit;
    q = h.predi;
    k--;
  } while (h.predi > 0);
}

BARRIER_INIT(barrier, NR_TASKLETS);
uint32_t packet_index = 0;
MUTEX_INIT(packet_index_lock);

void get_next_packet_index(uint32_t *index)
{
  mutex_lock(packet_index_lock);
  *index = packet_index;
  packet_index++;
  mutex_unlock(packet_index_lock);
}

void reset_packet_index()
{
  mutex_lock(packet_index_lock);
  packet_index = 0;
  mutex_unlock(packet_index_lock);
}

void decode(uint32_t packet_index_)
{

  /**
   * copy input sequence to WRAM buffer
   **/


  // packet_index_ = 0;
  uint64_t codetext_size = I.shapes[1];
  uint64_t codetext_aligned_size_byte = sizeof(Uchar) * I.aligned_shapes[1];

  __dma_aligned Uchar input_codetext__[MAX_CODETEXT_LEN];
  for (uint64_t i = 0; i < MAX_CODETEXT_LEN; i++)
    input_codetext__[i] = 0;
  mram_read(I.mram_addr[packet_index_], input_codetext__, codetext_aligned_size_byte);

#ifdef MESURE_BW
  if (perf_bw)
    nb_bytes_loaded[me()] += codetext_aligned_size_byte;
#endif

  reset_heap_and_stack();
  shoveltheheap((Int)(codetext_size), nmessbit, input_codetext__);
  // print_heap_final_pos();
  // printf("check DPU heap me () %u %d \n", me(), heap_check(0));
  


  /** reset vbits sequence **/
  uint64_t vbits_size;
  __dma_aligned Uchar vbits[MAX_DECODED_VBITS];
  for (uint64_t i = 0; i < MAX_DECODED_VBITS; i++)
    vbits[i] = 0;


  /** fill vbits vector from heap & stack of hypothesis **/
  traceback(vbits, &vbits_size);
  //return ;

  /** convert decoded vbits into regular message (decoded) bytes
   *  Both nmessbit (complete size of the decoded message in bytes),
   *  and vbits_size (effective size of the decoded vbits) to the function.
   *  vbits_size could be < to nmessbits, in this case part of the sequence
   *  won't be filled.
   **/
  __dma_aligned Uchar message[MAXSEQ_DPU];
  for (uint64_t i = 0; i < MAXSEQ_DPU; i++)
    message[i] = 0;
  packvbits(message, vbits, nmessbit, vbits_size);

  uint64_t vbits_aligned_size_byte = sizeof(Uchar) * O.aligned_shapes[1];
  mram_write(message, O.mram_addr[packet_index_], vbits_aligned_size_byte);

#ifdef MESURE_BW
  if (perf_bw)
    nb_bytes_written[me()] += vbits_aligned_size_byte;
#endif
}
// PROFILING_INIT(foo);

#if defined(ENABLE_HEAP_HOST_DEBUGGING) && ENABLE_HEAP_HOST_DEBUGGING == 1
uint64_t dbglimit = 200;
__dma_aligned heap_score_type dbgScores[200];
__dma_aligned heap_ptr_type dbgPtr[200];

#endif
int main()
{
  barrier_wait(&barrier);
  
  if (me() == 0)
  {
    assert((sizeof(hypothesis) % 8) == 0);

    if (perf_type == 2)
    {
      printf("CONFIG COUNT_INSTRUCTIONS \n");
      perfcounter_config(COUNT_INSTRUCTIONS, true);
    }
    else
    {
      printf("CONFIG COUNT_CYCLES \n");
      perfcounter_config(COUNT_CYCLES, true);
    }
  }
#if defined(MESURE_PERF)
  perfcounter_t start_time_;
  start_time_ = perfcounter_get();
#endif

  heap_init(&heap_, HEAP_MAX_ITEM, heap_items);
  barrier_wait(&barrier);

  if (me() == 0)
  {
    reward_ = DOUBLE_TO_FP(_reward);
    substitution_ = DOUBLE_TO_FP(_substitution);
    deletion_ = DOUBLE_TO_FP(_deletion);
    insertion_ = DOUBLE_TO_FP(_insertion);
    dither_ = DOUBLE_TO_FP(_dither);
    printf("reward_       %d\n", reward_);
    printf("substitution_ %d\n", substitution_);
    printf("deletion_     %d\n", deletion_);
    printf("insertion_    %d\n", insertion_);
    printf("dither_       %d\n", dither_);

    /**
     *  printf("HEAP MAX ITEM %u\n", (unsigned)(HEAP_MAX_ITEM));
     *  printf("TASKLETS %u\n", (unsigned)(NR_TASKLETS));
     *  printf("NSTAK %u\n", (unsigned)(NSTAK));
     *  printf("HLIMIT %u\n", (unsigned)(HLIMIT));
     *  printf("MAXSEQ_DPU %u\n", (unsigned)(MAXSEQ_DPU));
     *  printf("HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE %u\n", (unsigned)(HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE));
     *  printf("BUDDY_ALLOCATOR_SIZE_BYTE %u\n", (unsigned)(BUDDY_ALLOCATOR_SIZE_BYTE));
     **/
    xitf.init(&xitf, xferitf_buffer);

    if (!loaded)
    {
      // TODO : Not clear if buddy_init has to be done one time or not
      buddy_init(BUDDY_ALLOCATOR_SIZE_BYTE);

      /**
       * Get decoder patameters
       * **/
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

      /**
       * Get encoded input sequences
       * as Tensor2d.
       * The first dim defines the number
       * independant sequences (DNA strands) to decode.
       * The second dim defines the encoded sequence lenght.
       * **/
      xitf.getTensor2dUchar(&xitf, &I);

      /**
       * push output sequences as Tensor2d.
       * The first dim defines the number
       * independant decoded sequences.
       * The second dim defines the decoded sequence lenght.
       * NOTE : The basic elements of decoded sequences is bit,
       * but we use one byte (Uchar) to store each bit.
       * **/
      uint64_t nmessbytes = (nmessbit + 7) / 8;
      uint64_t output_shapes[2] = {I.shapes[0], nmessbytes};
      xitf.pushTensor2dUchar(&xitf, &O, output_shapes);

      codetextlen_g = I.shapes[1];

      /**
       * push heap (scores and ptr) tensor2d
       * if host debugging mode is activated
       * **/
      if (ENABLE_HEAP_HOST_DEBUGGING)
      {
        printf(" ENABLE_HEAP_HOST_DEBUGGING \n");
        uint64_t heap_shapes[2] = {1, HLIMIT};
        xitf.pushTensor2dUINT64(&xitf, &O_scores, heap_shapes);
        xitf.pushTensor2dUINT64(&xitf, &O_ptr, heap_shapes);
      }
      loaded = 1;
    }
    // return 0;

#if defined(MESURE_PERF)
    if (perf_type)
      nb_cycles_io[me()] += perfcounter_get() - start_time_;
#endif
    // check2dTensorUchar(&I);
    printf("perf type  %lu\n", perf_type);
  }

  printf("tasklet %u: stack = %u\n", me(), check_stack());
  printf("hyp  %u: \n", sizeof(hypothesis));
  printf("aligned hyp  %u: \n", XFER_MEM_ALIGN( sizeof(hypothesis)));
  barrier_wait(&barrier);
  uint32_t tasklet_current_packet_index;
  uint32_t total_packet = O.shapes[0];

  // profiling_start(&foo);

#if defined(MESURE_PERF)
  perfcounter_t start_time;
  if (perf_type)
    start_time = perfcounter_get();
#endif

  while (1)
  {
    get_next_packet_index(&tasklet_current_packet_index);
    if (tasklet_current_packet_index >= total_packet)
      break;
    decode(tasklet_current_packet_index);
  //  if(tasklet_current_packet_index >=10)
  //      break;
  }

#if defined(MESURE_PERF)
  if (perf_type)
    nb_cycles_decode[me()] += perfcounter_get() - start_time;
#endif

  //  profiling_stop(&foo);
  // barrier_wait(&barrier);

  if (perf_type)
  {
    total_cycles[me()] += perfcounter_get() - start_time_;
  }
#if defined(ENABLE_HEAP_HOST_DEBUGGING) && ENABLE_HEAP_HOST_DEBUGGING == 1
  // assert(NR_TASKLETS == 1);
  // mram_read(&heap_scores[me()][0], dbgScores, sizeof(Doub) * dbglimit);
  // mram_write(dbgScores, O_scores.mram_addr[me()], sizeof(Doub) * dbglimit);
  // mram_read(&heap_ptr[me()][0], dbgPtr, sizeof(uint64_t) * dbglimit);
  // mram_write(dbgPtr, O_ptr.mram_addr[me()], sizeof(uint64_t) * dbglimit);
#endif

  //#endif
  //
  //   if (me() == 0)
  //   {
  //     free2dTensor(&I);
  //     free2dTensor(&O);
  //     if (ENABLE_HEAP_HOST_DEBUGGING)
  //     {
  //       free2dTensor(&O_scores);
  //       free2dTensor(&O_ptr);
  //     }
  //   }

  return 0;
}

