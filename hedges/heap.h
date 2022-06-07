/**
 * @file heap.h
 * @author Dimitri Gerin (dgerin@upmem.com)
 * @brief HEDGES inner decoder scoring heap DPU implementation
 *        # original project #
 *        Based on original HEDGES project
 *        HEDGES Error-Correcting Code for DNA Storage Corrects Indels and
 *        Allows Sequence Constraints.
 *        William H. Press, John A. Hawkins, Stephen Knox Jones Jr,
 *        Jeffrey M. Schaub, Ilya J. Finkelstein
 *        submitted to Proceedings of the National Academy of Sciences.
 * @copyright 2022 UPMEM
 */

#ifndef __HEAP__
#define __HEAP__

#ifdef DPU_ENV
#include <assert.h>
#include <common.h>
#include <defs.h>
#include <mram.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#endif

typedef int32_t heap_score_type;
typedef uint32_t heap_ptr_type;

#ifdef DPU_ENV

/* perf counter */
#if defined(MESURE_PERF) || defined(MESURE_BW)
#include <perfcounter.h>
#endif
#include <profiling.h>
#ifdef MESURE_PERF
__host perfcounter_t nb_pop[NR_TASKLETS] = { 0 };
__host perfcounter_t nb_push[NR_TASKLETS] = { 0 };
__host perfcounter_t nb_cycles_pop[NR_TASKLETS] = { 0 };
__host perfcounter_t nb_cycles_push[NR_TASKLETS] = { 0 };
#endif
#ifdef MESURE_BW
__host uint64_t nb_bytes_loaded_heap[NR_TASKLETS] = { 0 };
__host uint64_t nb_bytes_written_heap[NR_TASKLETS] = { 0 };
#endif

/**
 * @brief Fixed Point Quantification macro  double -> fp
 */
#define HEAP_SCORE_MAXVAL_FLOAT 32767
#define FLOAT_TO_FP(i) (heap_score_type)((float)(i) * ((float)(1 << DECODER_QUANT_FRAC_BITS)))
#define FP_TO_FLOAT(i) (float)((float)(i) / ((float)(1 << DECODER_QUANT_FRAC_BITS)))

typedef struct item {
    heap_score_type score;
    heap_ptr_type ptr;
} item;

/**
 * @brief DPU heap datastructure
 *
 */
typedef struct heap {
    /* current position of the last heap element */
    __dma_aligned uint32_t heap_pos[NR_TASKLETS];
    /* items */
    __mram_ptr uint8_t *heap_items[NR_TASKLETS];
    /* max number of element */
    uint32_t size;
} heap;

/**
 * @brief initialize heap object
 *
 * @param h heap object pointer
 * @param size
 * @param heap_items prealocated array of item pointer
 * @return
 */
void heap_init(heap *h, uint32_t size_, __mram_ptr uint8_t *heap_items)
{

    if (!me()) {
        h->size = size_;
        assert(h->size > 0);

        for (uint8_t i = 0; i < NR_TASKLETS; i++) {
            h->heap_pos[i] = 0;
            h->heap_items[i] = heap_items + i * XFER_MEM_ALIGN(h->size * sizeof(item));
        }
    }
}

PROFILING_INIT(addr_compute);

/**
 * @brief read one item in heap
 *
 * @param h heap object pointer
 * @param h cur_pos position in the heap
 * @param it iterm pointer to fill
 * @return
 */
void heap_read(heap *h, uint32_t cur_pos, __dma_aligned item *it)
{

    __mram_ptr uint8_t *p = h->heap_items[me()] + cur_pos * sizeof(item);
    mram_read(p, it, sizeof(item));
#if defined(MESURE_BW) && (MESURE_BW == 1)
    nb_bytes_loaded_heap[me()] += sizeof(heap_ptr_type);
#endif
}
void heap_write(heap *h, uint32_t cur_pos, __dma_aligned item *it)
{
    __mram_ptr uint8_t *p = h->heap_items[me()] + cur_pos * sizeof(item);
    mram_write(it, p, sizeof(item));
#if defined(MESURE_BW) && (MESURE_BW == 1)
    nb_bytes_written_heap[me()] += sizeof(heap_ptr_type);
#endif
}

/**
 * @brief swap two items a and b in heap
 *
 * @param h heap object pointer
 * @param a position of item a
 * @param b position of item b
 * @return
 */
void heap_swap(heap *h, uint32_t a, uint32_t b)
{
    __dma_aligned item a_;
    __dma_aligned item b_;
    heap_read(h, a, &a_);
    heap_read(h, b, &b_);
    heap_write(h, a, &b_);
    heap_write(h, b, &a_);
}
/**
 * @brief heap push function
 *
 * @param h heap object pointer
 * @param it new item pointer
 * @param perf optional arg for perf measurement
 * @return
 */
void heap_push(heap *h, __dma_aligned item *it, bool perf)
{
#if defined(MESURE_BW) && (MESURE_BW == 1)
    nb_push[me()] += 1;
#endif

#if defined(MESURE_PERF)
    perfcounter_t start_time;
    if (perf)
        start_time = perfcounter_get();
#endif

    __dma_aligned item parent_it = (item) { .score = 0 };
    uint32_t parent_pos;
    uint32_t cur_pos;

    cur_pos = h->heap_pos[me()]++;

    /* store score and ptr at the last position */
    heap_write(h, cur_pos, it);

    assert(cur_pos < h->size);

    do {
        if (cur_pos <= 0)
            break;

        /* load parent node score */
        parent_pos = (cur_pos - 1) >> 1;
        heap_read(h, parent_pos, &parent_it);
        __dma_aligned item ro;
        heap_read(h, 0, &ro);
        if (parent_it.score > (*it).score) {
            /* swap parrent/child nodes */
            heap_swap(h, cur_pos, parent_pos);
            cur_pos = parent_pos;
        }
    } while (parent_it.score > (*it).score);

#if defined(MESURE_PERF)
    nb_cycles_push[me()] += perfcounter_get() - start_time;
#endif
}

/**
 * @brief heap pop function
 *
 * @param h heap object pointer
 * @param it item pointer to fill
 * @param perf optional arg for perf measurement
 * @return success
 */
bool heap_pop(heap *h, __dma_aligned item *it, bool perf)
{
#if defined(MESURE_BW) && (MESURE_BW == 1)
    nb_pop[me()] += 1;
#endif

#if defined(MESURE_PERF)
    perfcounter_t start_time;
    if (perf)
        start_time = perfcounter_get();
#endif

    uint32_t end_pos;
    uint32_t cur_pos;
    __dma_aligned item cur_item;

    uint32_t right_pos;
    uint32_t left_pos;
    uint32_t min_pos;
    __dma_aligned item left_item;
    __dma_aligned item right_item;
    heap_score_type min_score;

    __dma_aligned item maxitem = (item) { .score = FLOAT_TO_FP(HEAP_SCORE_MAXVAL_FLOAT), 0 };

    cur_pos = 0;

    if (h->heap_pos[me()] == 0)
        return false;

    /* assign score (poped score) at the top of the heap */
    heap_read(h, 0, it);

    /* init end_pos at las heap decremented position */
    end_pos = --(h->heap_pos[me()]);

    if (end_pos > 0) {
        /*
         * swap end and top nodes
         */
        heap_swap(h, 0, end_pos);

        /*
         * assign the end nodes (score and ptr) to specific values
         */
        heap_write(h, end_pos, &maxitem);
        /* NOTE : Not sure this step is necessary (only asign score value to
         * bigvalue should be necessary) heap_write_ptr(end_pos, &bigval);
         */

        /*
         * visit heap from top to last heap position (end_pos)
         * start by left child
         * compare left and right child to select the visiting
         * direction (left or right)
         */
        while ((left_pos = (cur_pos << 1) + 1) < end_pos) {
            /* read cur_score */
            heap_read(h, cur_pos, &cur_item);
            /* read left and right child */
            right_pos = left_pos + 1;
            heap_read(h, left_pos, &left_item);
            heap_read(h, right_pos, &right_item);

            /* find the child node with the minimal score */
            min_pos = (left_item.score < right_item.score ? left_pos : right_pos);
            min_score = (left_item.score < right_item.score ? left_item.score : right_item.score);

            if (cur_item.score > min_score) {
                heap_swap(h, cur_pos, min_pos);
            } else
                break;
            cur_pos = min_pos;
        }
    }

#if defined(MESURE_PERF)
    if (perf)
        nb_cycles_pop[me()] += perfcounter_get() - start_time;
#endif

    return true;
}

#endif

#endif
