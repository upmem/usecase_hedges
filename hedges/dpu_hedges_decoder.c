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


/**
 * @brief section performances counters
 */
#if defined(MESURE_PERF)
__host perfcounter_t nb_cycles_total = 0;
__host perfcounter_t nb_cycles_hypcompute[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_decode[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_io[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_hashfunc[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_penality[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_hypcompute_first_section[NR_TASKLETS] = {0};
__host perfcounter_t nb_cycles_hypload[NR_TASKLETS] = {0};
__host perfcounter_t total_cycles[NR_TASKLETS] = {0};
#endif

/**
 * @brief bandwidth counters
 */
#ifdef MESURE_BW
__host uint64_t nb_bytes_loaded[NR_TASKLETS] = {0};
__host uint64_t nb_bytes_written[NR_TASKLETS] = {0};
#endif

/**
 * @brief representation of HEDGES decoding hypothesis 
 * 
 */
typedef struct hypothesis{
  /** next char in message */
  Int offset;
  /** my position in the decoded message (0,1,...) */
  Int seq;
  /** last decoded up to now */
  Mbit messagebit;
  /** corresponding salt values */
  Ullong prevbits, salt, newsalt;
  /** precedent code */
  GF4reg prevcode;
  /** cumulated score */
  heap_score_type score;
  /** index of predecessor in hypostackp */
  Int predi;
} hypothesis;

__dma_aligned GF4reg acgtacgt;

/**
 * @brief initialize hypothesis node
 * 
 * @param h 
 * @return * cache* 
 */
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

/**
 * @brief tasklet wise stack of hypothesis
 */
__mram_noinit uint8_t hypothesis_stack[NR_TASKLETS][XFER_MEM_ALIGN(sizeof(hypothesis) * HEAP_MAX_ITEM)];

/**
 * @brief tasklet wise heao of referenced hypothesis
 */
__dma_aligned heap hypothesis_heap;
__mram_noinit uint8_t hypothesis_heap_buffer[NR_TASKLETS * XFER_MEM_ALIGN(sizeof(item) * HEAP_MAX_ITEM)];
__dma_aligned uint64_t heap_pos[NR_TASKLETS];
__dma_aligned uint64_t heap_pos_final[NR_TASKLETS];
__dma_aligned uint64_t nhypo[NR_TASKLETS];

/** first DPU run reference counter **/
__host uint64_t first_run = 1;
/** persistant xferIft offset snapshot **/
__host uint64_t xferOffset;
/** host accessible NR_TASKLETS value **/
__host uint64_t nr_tasklets = NR_TASKLETS;


/**
 * @brief HOST/DPU xferItf MRAM buffer
 */
__mram_noinit uint8_t xferitf_buffer[1 << HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE];


__host uint64_t perf_type;
__host uint64_t perf_bw;
__host uint64_t nmessbit;

/** constants seted at programm start */
Int codetextlen_g;

/**
 * @brief HOST PROVIDED (only for _reward) penality values
 * defined in double and quantified at programm first_run
 */
__host double _reward = -0.13;
double _substitution = 1.;
double _deletion = 1.;
double _insertion = 1.;
double _dither = 0.;
heap_score_type reward_;
heap_score_type substitution_;
heap_score_type deletion_;
heap_score_type insertion_;
heap_score_type dither_;

/**
 * @brief Fixed Point Quantification macro  double -> fp
 */
#define DOUBLE_TO_FP(i) (heap_score_type)((double)(i) * ((double)(1 << DECODER_QUANT_FRAC_BITS)))

/**
 * @brief HEDGES decoder parameters
 */

/** Greedy Exhaustive Search Maximal number of hypothesis for each strand to decode */
__host Int HLIMIT;
__host Int NSTAK;

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
__dma_aligned uint64_t RIGHTPRIMER_LEN;
__dma_aligned uint64_t PRIMERSALT_LEN;
__dma_aligned uint64_t PATTARR_LEN;
__dma_aligned GF4word rightprimer;
__dma_aligned VecUllong primersalt;
__dma_aligned VecInt pattarr;

/**
 *  IO Tensors : Tensor class provides  API for HOST/DPU exchanges
 *               of tensors. 
 *               I provides a batch of strands to decode from the HOST
 *               O provides the corresponding decoded sequence of bits
 */
Tensor2d I;
Tensor2d O;

#if (ENABLE_HEAP_HOST_DEBUGGING == 1)
Tensor2d O_scores, O_ptr;
#endif

/**
 * @brief returns the number of allowed ACGTs and puts them in dnac_ok
 * 
 * @param dnac_ok 
 * @param prev 
 * @return returns 
 */
Int dnacallowed(Uchar *dnac_ok, GF4reg prev)
{
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

/**
 * @brief encoder/decoder hash function 
 * 
 * @param out 
 * @param bits 
 * @param seq 
 * @param salt 
 * @param mod 
 */
void hash_dna(Int *out, Ullong bits, Int seq, Ullong salt, Int mod)
{
#if defined(MESURE_PERF)
  perfcounter_t start_time_______;
  if (perf_type)
    start_time_______ = perfcounter_get();
#endif

  Ullong ret;

  ranhash_int64(&ret, ((((((Ullong)(seq)&seqnomask) << NPREV) | bits) << HSALT) | salt));

  Ullong div;
  Ullong ret_;


  /* op: (% mod) with mod from 0 to 4 */
  if (mod == 4)
  {
    div = ret >> 2;
    ret_ = ret - (div << 2);
  }
  else if(mod == 2)
  {
    div = ret >> 1;
    ret_ = ret - (div << 1);
  }
  else if (mod == 1)
  {
    ret_ = ret;
  }
  else
  {
    ret_ = ret % mod;
  }

  *out = (Int)(ret_);
#if defined(MESURE_PERF)
 if (perf_type)
   nb_cycles_hashfunc[me()] += perfcounter_get() - start_time_______ ;
#endif
}

/**
 * @brief  Computes the message len from vbits mapping patern.
 * @param nmb 
 * @return Int 
 */
Int vbitlen(Int nmb)
{
  Int ksize, nn = 0;
  for (ksize = 0;; ksize++)
  {
    if (nn >= nmb)
      break;
    assert(ksize < ENCODER_STACK_MAXSEQ_DPU_NBYTES && "vbitlen: MAXSEQ too small");
    nn += pattarr[ksize];
  }
  return ksize;
}

/**
 * @brief transform decoded bit sequence info efficient packed bit sequence stored 
 * 
 * @param packed_bytes 
 * @param vbits 
 * @param nmessbits 
 * @param vbits_size 
 */
void packvbits(Uchar *packed_bytes, Uchar *vbits, uint64_t nmessbits, uint64_t vbits_size)
{
  Int i, j, k, k1;
  Uchar bit;
  uint64_t nn;
  assert(vbits_size <= MAX_DECODED_VBITS);

  /** number of bits */
  for (k = 0; k < vbits_size; k++)
    nn += pattarr[k];
  /** no more than the specified number of bits */
  nn = MIN(nn, nmessbits);
  nn = nmessbits;
  /** number of bytes */
  uint64_t nbytes = (nn + 7) / 8;

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



void hypothesis_load(__dma_aligned hypothesis *cache,  uint32_t pos)
{
  mram_read(&(hypothesis_stack[me()][pos * sizeof(hypothesis) ]), cache, XFER_MEM_ALIGN(sizeof(hypothesis)));
#ifdef MESURE_BW
  if (perf_bw)
    nb_bytes_loaded[me()] += sizeof(hypothesis);
#endif
}

void hypothesis_store(__dma_aligned hypothesis *cache, uint32_t pos)
{
  mram_write(cache, &(hypothesis_stack[me()][pos * sizeof(hypothesis)]), XFER_MEM_ALIGN(sizeof(hypothesis)));
#ifdef MESURE_BW
  if (perf_bw)
    nb_bytes_written[me()] += sizeof(hypothesis);
#endif
}

/**
 * @brief global reset of all hypothesis structures 
 */
void reset_hypothesis()
{
  __dma_aligned heap_ptr_type zero = 0;
  __dma_aligned hypothesis cached_h;
  __dma_aligned item cur_heap_item;

  hypothesis_init_root(&cached_h);
  hypothesis_store(&cached_h,  0);

  /** set nhypo to 1 element */
  nhypo[me()] = 1;
  hypothesis_heap.heap_pos[me()] = 0;

  cur_heap_item.score = 1000000000;
  cur_heap_item.ptr = 0;
  heap_push(&hypothesis_heap, &cur_heap_item, perf_type);
  return;
}

/**
 * @brief create new hypothesis from parent hyporhesis
 * 
 * @param h 
 * @param pred 
 * @param mbit 
 * @param skew 
 * @param input_codetext 
 * @return Int 
 */
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

  /** load precedent hypothesis */
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
  h->messagebit = mbit; /** variable number */
  h->seq = hp.seq + 1;
  if (h->seq > ENCODER_STACK_MAXSEQ_DPU_NBYTES)
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
    h->newsalt = ((hp.newsalt << 1) & saltmask) ^ h->messagebit; /** variable bits overlap, but that's ok with XOR */
  }
  else if (h->seq == NSP)
  {
    /** update salt value */
    mysalt = h->salt = hp.newsalt;
  }
  else
    mysalt = h->salt;
  h->offset = hp.offset + 1 + skew;
  if (h->offset >= codetextlen_g)
    return 0;

  /** calculate predicted message under this hypothesis */
  h->prevcode = hp.prevcode;
  mod = (h->seq < LPRIMER ? 4 : dnacallowed(dnac_ok, h->prevcode));

#if defined(MESURE_PERF)
  if (perf_type)
    nb_cycles_hypcompute_first_section[me()] += perfcounter_get() - start_time___;
#endif

   hash_dna(&regout, h->prevbits, h->seq, mysalt, mod);

#if defined(MESURE_PERF)
  perfcounter_t start_time_____;
  if (perf_type)
    start_time_____  = perfcounter_get();
#endif
  regout = (regout + (Uchar)(h->messagebit)) % mod;
  regout = (h->seq < LPRIMER ? regout : dnac_ok[regout]);
  h->prevbits = ((hp.prevbits << nbits) & prevmask) | h->messagebit; // variable number
  h->prevcode = ((h->prevcode << 2) | regout) & dnawinmask;


  /** deletion */
  if (skew < 0)
  {
    mypenalty = deletion_;
  }
  else
  {
    /** the only place where a check is possible!  **/
    discrep = (regout == input_codetext[h->offset]);
    /** substitution of error free */
    if (skew == 0)
      mypenalty = (discrep ? reward_ : substitution_);
    else
    {
      /** inserion */
      mypenalty = insertion_ + (discrep ? reward_ : substitution_);
    }
  }

#if defined(MESURE_PERF)
  if (perf_type)
    nb_cycles_penality[me()] += perfcounter_get() - start_time_____;
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

/**
 * @brief Principal decoder function that perform HED(GES) Greedy Exhaustive Search algorithm
 * 
 * @param limit 
 * @param nmessbits 
 * @param input_codetext 
 * @return
 */
void heap_mining(Int limit, Int nmessbits, __dma_aligned Uchar *input_codetext)
{
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
    bool nemptyheap = heap_pop(&hypothesis_heap, &cur_heap_item, perf_type);

    assert(nemptyheap && "pop empty heap");

    /** load hypothesis from stack */
    hypothesis_load(&hp,  cur_heap_item.ptr);

    seq = hp.seq;
    assert(seq < ENCODER_STACK_MAXSEQ_DPU_NBYTES && "heap_mining: MAXSEQ too small");

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

    /** sequential computing of each kind of Hypothesis, for each bits */
    /** error free hypothesis */
    for (mbit = 0; mbit < nguess; mbit++)
    {
      if (init_from_predecessor(&new_hp, cur_heap_item.ptr , mbit, 0, input_codetext))
      {
        new_heap_item.score = new_hp.score;
        new_heap_item.ptr = nhypo[me()];

        heap_push(&hypothesis_heap, &new_heap_item, perf_type);
        hypothesis_store(&new_hp , nhypo[me()]);
        nhypo[me()]++;
      }
    }
    /** deletion error hypothesis */
    for (mbit = 0; mbit < nguess; mbit++)
    {
      if (init_from_predecessor(&new_hp, cur_heap_item.ptr , mbit, -1, input_codetext))
      {
        new_heap_item.score = new_hp.score;
        new_heap_item.ptr = nhypo[me()];

        heap_push(&hypothesis_heap, &new_heap_item, perf_type);
        hypothesis_store(&new_hp, nhypo[me()]);
        nhypo[me()]++;
      }
    }
    /** insertion error hypothesis */
    for (mbit = 0; mbit < nguess; mbit++)
    {
      if (init_from_predecessor(&new_hp, cur_heap_item.ptr , mbit, 1, input_codetext))
      {
        new_heap_item.score = new_hp.score;
        new_heap_item.ptr = nhypo[me()];

        heap_push(&hypothesis_heap, &new_heap_item, perf_type);
        hypothesis_store(&new_hp, nhypo[me()]);
        nhypo[me()]++;
      }
    }
  }
  heap_pos_final[me()] = cur_heap_item.ptr;
}

/**
 * @brief trace the heap hypothesis back to extract the decoded bit sequence 
 * 
 * @param limit 
 * @param nmessbits 
 * @param input_codetext 
 * @return
 */
void heap_traceback(__dma_aligned Uchar *vbits, uint64_t *traceback_size)
{

  Int k, kk = 0;
  Int nfinal = (Int)(heap_pos_final[me()]);
  Int q = nfinal;
  __dma_aligned hypothesis h;
  /**
   * Computes the length of the decoded bits sequence
   * starting from the heap final position.
   * **/
  do
  {
    hypothesis_load(&h,  q);
    if (!(h.predi > 0))
      break;
    kk++;
    q = h.predi;
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

/** principal tasklet barrier */
BARRIER_INIT(barrier, NR_TASKLETS);
/** representation of curently processed packet index (TASKLET shared). */
uint32_t packet_index = 0;
/** tasklet lock for packet_index RW **/
MUTEX_INIT(packet_index_lock);

void get_next_packet_index(uint32_t *index)
{
  mutex_lock(packet_index_lock);
  *index = packet_index;
  packet_index++;
  mutex_unlock(packet_index_lock);
}

/**
 * @brief reset given packet index 
 * 
 * @param limit 
 * @param nmessbits 
 * @param input_codetext 
 * @return
 */
void reset_strands_packet_index()
{
  mutex_lock(packet_index_lock);
  packet_index = 0;
  mutex_unlock(packet_index_lock);
}

/**
 * @brief decode one input sequence (from ACGT space to binary space)
 * 
 * @param limit 
 * @param nmessbits 
 * @param input_codetext 
 * @return
 */
void decode(uint32_t packet_index_)
{

  /**
   * copy input sequence to WRAM buffer
   **/
  uint64_t codetext_size = I.shapes[1];
  uint64_t codetext_aligned_size_byte = sizeof(Uchar) * I.aligned_shapes[1];

  __dma_aligned Uchar input_codetext__[ENCODER_STACK_MAXSEQ_DPU_NBYTES];
  for (uint64_t i = 0; i < ENCODER_STACK_MAXSEQ_DPU_NBYTES; i++)
    input_codetext__[i] = 0;
  mram_read(I.mram_addr[packet_index_], input_codetext__, codetext_aligned_size_byte);

// #if (MESURE_BW ==1)
//   if (perf_bw)
//     nb_bytes_loaded[me()] += codetext_aligned_size_byte;
// #endif

  reset_hypothesis();
  heap_mining((Int)(codetext_size), nmessbit, input_codetext__);
  
#if (ENABLE_HEAP_HOST_DEBUGGING == 1)
  printf("[HEDGES][DECODER][HEAP FINAL POS] tasklet[%u], packet [%lu] : %lu \n", (unsigned)me(), packet_index_, heap_pos_final[me()]);
#endif

  /** reset vbits sequence **/
  uint64_t vbits_size;
  __dma_aligned Uchar vbits[MAX_DECODED_VBITS];
  for (uint64_t i = 0; i < MAX_DECODED_VBITS; i++)
    vbits[i] = 0;


  /** fill vbits vector from heap & stack of hypothesis **/
  heap_traceback(vbits, &vbits_size);

  /** convert decoded vbits into regular message (decoded) bytes
   *  Both nmessbit (complete size of the decoded message in bytes),
   *  and vbits_size (effective size of the decoded vbits) to the function.
   *  vbits_size could be < to nmessbits, in this case part of the sequence
   *  won't be filled.
   **/
  __dma_aligned Uchar message[ENCODER_STACK_MAXSEQ_DPU_NBYTES];
  for (uint64_t i = 0; i < ENCODER_STACK_MAXSEQ_DPU_NBYTES; i++)
    message[i] = 0;
  packvbits(message, vbits, nmessbit, vbits_size);

  uint64_t vbits_aligned_size_byte = sizeof(Uchar) * O.aligned_shapes[1];
  mram_write(message, O.mram_addr[packet_index_], vbits_aligned_size_byte);

#if (MESURE_BW == 1)
  if (perf_bw)
    nb_bytes_written[me()] += vbits_aligned_size_byte;
#endif
}

#if (ENABLE_HEAP_HOST_DEBUGGING == 1)
uint64_t dbglimit = 200;
__dma_aligned heap_score_type dbgScores[200];
__dma_aligned heap_ptr_type dbgPtr[200];
#endif

/**
 * @brief DPU HEDGES decoder entry point 
 * 
 * @param limit 
 * @param nmessbits 
 * @param input_codetext 
 * @return
 */
int main()
{


  barrier_wait(&barrier);
  
  if (me() == 0)
  {
    assert((sizeof(hypothesis) % 8) == 0);

    if (perf_type == 2)
    {
      printf("[HEDGES][DECODER][PERFCOUNTER] INSTRUCTIONS\n");
      perfcounter_config(COUNT_INSTRUCTIONS, true);
    }
    else
    {
      printf("[HEDGES][DECODER][PERFCOUNTER] CYCLES\n");
      perfcounter_config(COUNT_CYCLES, true);
    }
  }
#if (MESURE_PERF == 1)
  perfcounter_t start_time_;
  start_time_ = perfcounter_get();
#endif

  heap_init(&hypothesis_heap, HEAP_MAX_ITEM, hypothesis_heap_buffer);
  barrier_wait(&barrier);

  if (me() == 0)
  {
    reward_ = DOUBLE_TO_FP(_reward);
    substitution_ = DOUBLE_TO_FP(_substitution);
    deletion_ = DOUBLE_TO_FP(_deletion);
    insertion_ = DOUBLE_TO_FP(_insertion);
    dither_ = DOUBLE_TO_FP(_dither);

    printf("[HEDGES][DECODER][QUANT] reward        %d\n", reward_);
    printf("[HEDGES][DECODER][QUANT] substitution  %d\n", substitution_);
    printf("[HEDGES][DECODER][QUANT] deletion      %d\n", deletion_);
    printf("[HEDGES][DECODER][QUANT] insertion     %d\n", insertion_);
    printf("[HEDGES][DECODER][QUANT] dither        %d\n", dither_);

    xitf.init(&xitf, xferitf_buffer);

    if (first_run)
    {
      /**
       * Global init of buddy allocator (used by xferItf)
       * **/
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
      /**
       * save a snapshot of xferItf buffer Offset
       * **/
      xitf.save(&xitf, &xferOffset);

      /**
       * restore the snapshot of xferItf buffer Offset
       * **/
      xitf.restore(&xitf, &xferOffset);

      /**
       * reset strands Packet Index
       * **/
      reset_strands_packet_index();

      /**
       * get encoded input sequences
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
#if (ENABLE_HEAP_HOST_DEBUGGING == 1)
        printf(" PUSH HEAP TO HOST\n");
        uint64_t heap_shapes[2] = {1, HLIMIT};
        xitf.pushTensor2dUINT64(&xitf, &O_scores, heap_shapes);
        xitf.pushTensor2dUINT64(&xitf, &O_ptr, heap_shapes);
#endif

      first_run = 1;
    }

#if (MESURE_PERF == 1)
    if (perf_type)
      nb_cycles_io[me()] += perfcounter_get() - start_time_;
#endif
  }
  printf("[HEDGES][DECODER][INFO] tasklet %u: stack_size  %u Byte\n", me(), check_stack());

  barrier_wait(&barrier);
  uint32_t tasklet_current_packet_index;
  uint32_t total_packet = O.shapes[0];

#if (MESURE_PERF == 1)
  perfcounter_t start_time;
  if (perf_type)
    start_time = perfcounter_get();
#endif

  /**
   * @brief decode packet per packet
   **/
  do
  {
    get_next_packet_index(&tasklet_current_packet_index);
    if (tasklet_current_packet_index >= total_packet)
      break;
    decode(tasklet_current_packet_index);
  } while (1);

#if (MESURE_PERF == 1)
  if (perf_type)
    nb_cycles_decode[me()] += perfcounter_get() - start_time;
#endif

  if (perf_type)
  {
    total_cycles[me()] += perfcounter_get() - start_time_;
  }


  nb_cycles_total = perfcounter_get() ;

  printf("nb_cycles_total() %lu\n", perfcounter_get() );

  return 0;
}

