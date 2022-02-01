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
#include <profiling.h>
#include "heap.h"

__mram_noinit uint8_t xferitf_buffer[1 << HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE];
__host __dma_aligned uint64_t loaded = 0;
__host __dma_aligned uint64_t xferOffset;

// MAXSEQ_DPU is used to set buffer size of both unpacked
// bits buffer, and uncoded bit buffer. Theses two buffer
// are declared in TASKLET stack, so we must limit its size
// NOTE that the effective size of both buffer are computed as
// uint64_t codetext_size = vbitlen(8 /*8 bit per packed bit*/ * packed_message_size) + RPRIMER;
// but we must use constants for C array SIZE, so we fix the limit
// at the smallest possible value right now.
// WARNING : 300 works only for coderate = 0.5
#define MAXSEQ_DPU ENCODER_STACK_MAXSEQ_DPU_NBYTES

__host uint64_t perf_count;
__host uint64_t perf_cycles_or_inst;
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
__host double reward = -0.13;
double substitution = 1.;
double deletion = 1.;
double insertion = 1.;
double dither = 0.;

void quantif_test()
{
  double a = 0.0;
  double b;
  b = a;
  for (uint64_t i = 0; i < 10; i++)
  {
    b = b + reward;
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

// PROFILING_INIT(foo);

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
  return (Int)(ranhash_int64((((((Ullong)(seq)&seqnomask) << NPREV) | bits) << HSALT) | salt) % mod);
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
  printf(" nbytes %lu, vbits size %lu \n", nbytes, vbits_size);

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
  __dma_aligned double score;
  /** index of predecessor in hypostackp **/
  __dma_aligned Int predi;
} hypothesis;

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

__mram_noinit hypothesis hstack[NR_TASKLETS][HEAP_MAX_ITEM];

void hypothesis_load(__dma_aligned hypothesis *cache, uint16_t tasklet_id, uint32_t pos)
{
  mram_read(&(hstack[tasklet_id][pos]), cache, sizeof(hypothesis));
}

void print_hypothesis(uint16_t tasklet_id, uint32_t pos)
{
  hypothesis cached_h;
  hypothesis_load(&cached_h, me(), pos);
  printf(" predi %d  \n", cached_h.predi);
  printf(" offset %d  \n", cached_h.offset);
  printf(" seq %d  \n", cached_h.seq);
  printf(" messagebit %d  \n", cached_h.messagebit);
  printf(" prevbits %llu  \n", cached_h.prevbits);
  printf(" score %lf  \n", cached_h.score);
  printf(" salt %llu  \n", cached_h.salt);
  printf(" newsalt %llu  \n", cached_h.newsalt);
  printf(" prevcode %llu  \n", cached_h.prevcode);
}

void hypothesis_store(__dma_aligned hypothesis *cache, uint16_t tasklet_id, uint32_t pos)
{
  mram_write(cache, &(hstack[tasklet_id][pos]), sizeof(hypothesis));
}

void init_heap_and_stack()
{
  __dma_aligned hypothesis cached_h;

  heap_init();

  hypothesis_init_root(&cached_h);
  hypothesis_store(&cached_h, me(), 0);

  /** set nhypo to 1 element **/
  nhypo[me()] = 1;

  heap_pos[me()] = 0;
  // heap_pos_final[me()] = 0;

  __dma_aligned double score = 10000000000;
  __dma_aligned uint64_t ptr = 0;
  heap_push(&score, &ptr);
  return;
}

Int init_from_predecessor(hypothesis *h, Int pred, Mbit mbit, Int skew, __dma_aligned Uchar *input_codetext)
{
  bool discrep;
  Int regout, mod;
  double mypenalty;
  Ullong mysalt;

  /** load precedent hypothesis**/
  hypothesis hp;
  hypothesis_load(&hp, me(), pred);

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
  regout = digest(h->prevbits, h->seq, mysalt, mod);
  regout = (regout + (Uchar)(h->messagebit)) % mod;
  regout = (h->seq < LPRIMER ? regout : dnac_ok[regout]);
  h->prevbits = ((hp.prevbits << nbits) & prevmask) | h->messagebit; // variable number
  h->prevcode = ((h->prevcode << 2) | regout) & dnawinmask;

  /**deletion**/
  if (skew < 0)
  {
    mypenalty = deletion;
  }
  else
  {
    /** the only place where a check is possible!  **/
    discrep = (regout == input_codetext[h->offset]);
    /** substitution of error free**/
    if (skew == 0)
      mypenalty = (discrep ? reward : substitution);
    else
    {
      /** inserion **/
      mypenalty = insertion + (discrep ? reward : substitution);
    }
  }

  /** NOTE : dont't understant why this code does.
   * It adds a randomly signed zero -0./+O.
   * TOTO : understand the impact of this step
   * in the scoring.
   *
   * if (dither > 0.)
   *  mypenalty += dither * (2. * ran_double() - 1.);
   **/
  h->score = hp.score + mypenalty;
  return 1;
}

void shoveltheheap(Int limit, Int nmessbits, __dma_aligned Uchar *input_codetext)
{
  // given the heap, keep processing it until offset limit, hypothesis limit, or an error is reached
  Int seq, nguess, qqmax = -1, ofmax = -1, seqmax = vbitlen(nmessbits);
  Uchar mbit;
  __dma_aligned double curscore;
  __dma_aligned uint64_t qq;
  __dma_aligned hypothesis hp;
  uint8_t errcode = 0;
  while (true)
  {
    /**
     * pop heap -> gives best curscore and its associated
     * hypothesis position (hypothesis stack)
     * **/
    bool nemptyheap = heap_pop(&curscore, &qq);

    assert(nemptyheap && "pop empty heap");

    /**load hypothesis from stack**/
    hypothesis_load(&hp, me(), qq);

    seq = hp.seq;
    assert(seq < MAXSEQ_DPU && "shoveltheheap: MAXSEQ too small");

    nguess = 1 << pattarr[seq + 1]; // i.e., 1, 2, or 4
    if (hp.offset > ofmax)
    { // keep track of farthest gotten to
      ofmax = hp.offset;
      qqmax = qq;
    }
    if (curscore > 1.e10)
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

    /**
     * sequential computing of each kind of Hypothesis, for each bits
     * **/
    /** error free hypothesis **/
    for (mbit = 0; mbit < nguess; mbit++)
    {
      hypothesis new_hp;
      if (init_from_predecessor(&new_hp, qq, mbit, 0., input_codetext))
      {
        heap_push(&(new_hp.score), &nhypo[me()]);
        hypothesis_store(&new_hp, me(), nhypo[me()]);
        nhypo[me()]++;
      }
    }
    /** deletion error hypothesis **/
    for (mbit = 0; mbit < nguess; mbit++)
    {
      hypothesis new_hp;
      if (init_from_predecessor(&new_hp, qq, mbit, -1., input_codetext))
      {
        heap_push(&(new_hp.score), &nhypo[me()]);
        hypothesis_store(&new_hp, me(), nhypo[me()]);
        nhypo[me()]++;
      }
    }
    /** insertion error hypothesis **/
    for (mbit = 0; mbit < nguess; mbit++)
    {
      hypothesis new_hp;
      if (init_from_predecessor(&new_hp, qq, mbit, 1., input_codetext))
      {
        heap_push(&(new_hp.score), &nhypo[me()]);
        hypothesis_store(&new_hp, me(), nhypo[me()]);
        nhypo[me()]++;
      }
    }
  }
  heap_pos_final[me()] = qq;
}

void print_heap_final_pos()
{
  printf("================> heap_final_pos %lu \n", heap_pos_final[me()]);
}

/**
 * extract vbits from heap
 * **/
void traceback(__dma_aligned Uchar *vbits, uint64_t *traceback_size)
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
    hypothesis_load(&h, me(), q);
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
    hypothesis_load(&h, me(), q);
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
  init_heap_and_stack();
  shoveltheheap((Int)(codetext_size), nmessbit, input_codetext__);
  // printf("check DPU heap me () %u %d \n", me(), heap_check(0));

  print_heap_final_pos();

  /** reset vbits sequence **/
  uint64_t vbits_size;
  __dma_aligned Uchar vbits[MAX_DECODED_VBITS];
  for (uint64_t i = 0; i < MAX_DECODED_VBITS; i++)
    vbits[i] = 0;

  /** fill vbits vector from heap & stack of hypothesis **/
  traceback(vbits, &vbits_size);

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
}

int main()
{

  if (me() == 0)
  {
    if (perf_cycles_or_inst)
    {
      perfcounter_config(COUNT_CYCLES, true);
      printf(" count CYCLES ?  \n");
    }
    else
    {
      perfcounter_config(COUNT_INSTRUCTIONS, true);
      printf(" count INST ?  \n");
    }
    printf("HEAP MAX ITEM %u\n", (unsigned)(HEAP_MAX_ITEM));
    printf("TASKLETS %u\n", (unsigned)(NR_TASKLETS));
    printf("NSTAK %u\n", (unsigned)(NSTAK));
    printf("HLIMIT %u\n", (unsigned)(HLIMIT));
    printf("MAXSEQ_DPU %u\n", (unsigned)(MAXSEQ_DPU));
    printf("HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE %u\n", (unsigned)(HOSTDPU_XFER_BUFFER_SIZE_LOG2_BYTE));
    printf("BUDDY_ALLOCATOR_SIZE_BYTE %u\n", (unsigned)(BUDDY_ALLOCATOR_SIZE_BYTE));
    assert(64 == sizeof(hypothesis));
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
        uint64_t heap_shapes[2] = {1, HLIMIT};
        xitf.pushTensor2dUINT64(&xitf, &O_scores, heap_shapes);
        xitf.pushTensor2dUINT64(&xitf, &O_ptr, heap_shapes);
      }
      loaded = 1;
    }
    // check2dTensorUchar(&I);
  }

  printf("tasklet %u: stack = %u\n", me(), check_stack());
  barrier_wait(&barrier);
  uint32_t tasklet_current_packet_index;
  uint32_t total_packet = O.shapes[0];

  while (1)
  {
    get_next_packet_index(&tasklet_current_packet_index);
    if (tasklet_current_packet_index >= total_packet)
      break;
    decode(tasklet_current_packet_index);
  }

  barrier_wait(&barrier);
  if (ENABLE_HEAP_HOST_DEBUGGING)
  {
    /** copy heap to xferItf correspondant Tensor2d **/
    for (uint64_t i = 0; i < HLIMIT; i++)
    {
      double sc;
      __mram_ptr double *score = ((__mram_ptr double *)(heap_scores[me()]) + i);
      __mram_ptr double *score_b = ((__mram_ptr double *)(O_scores.mram_addr[0]) + i);
      mram_read(score, &sc, 8);
      mram_write(&sc, score_b, 8);
      uint64_t ptrr;
      __mram_ptr double *ptr = ((__mram_ptr double *)(heap_ptr[me()]) + i);
      __mram_ptr uint64_t *ptr_b = ((__mram_ptr uint64_t *)(O_ptr.mram_addr[0]) + i);
      mram_read(ptr, &ptrr, 8);
      mram_write(&ptrr, ptr_b, 8);
    }
  }

  if (me() == 0)
  {
    free2dTensor(&I);
    free2dTensor(&O);
    if (ENABLE_HEAP_HOST_DEBUGGING)
    {
      free2dTensor(&O_scores);
      free2dTensor(&O_ptr);
    }
    perf_count = perfcounter_get();
  }

  printf("DONE\n");
  return 0;
}
