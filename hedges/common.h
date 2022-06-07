/*
 * @file common.h
 * @author Dimitri Gerin (dgerin@upmem.com)
 * @brief host/dpu shared common Header file (CPU side)
 *        # original project #
 *        Based on original HEDGES project
 *        HEDGES Error-Correcting Code for DNA Storage Corrects Indels and
 *        Allows Sequence Constraints.
 *        William H. Press, John A. Hawkins, Stephen Knox Jones Jr,
 *        Jeffrey M. Schaub, Ilya J. Finkelstein
 *        submitted to Proceedings of the National Academy of Sciences.
 * @copyright 2022 UPMEM
 */

#ifndef __COMMON__
#define __COMMON__
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#ifndef DPU_ENV
#include <vector>
#else
#include <alloc.h>
#include <mram.h>
#endif

/* 8 byte allignement helpers */
#define MEMORY_ALIGNMENT_BYTE_LOG2 3
#define MEMORY_ALIGNMENT_BYTE (1 << MEMORY_ALIGNMENT_BYTE_LOG2)
#define XFER_MEM_ALIGN(expr)                                                                                                     \
    (((((expr) + MEMORY_ALIGNMENT_BYTE - 1) >> MEMORY_ALIGNMENT_BYTE_LOG2)) << MEMORY_ALIGNMENT_BYTE_LOG2)

/* size of the buffer for which we compute the checksum 64 KB */
#define BUFFER_SIZE (1 << BUFFER_SIZE_LOG2)

/* basic type names (redefine if your bit lengths don't match) */
/* 32 bit integer */
typedef int Int;
typedef unsigned int Uint;
/* 8 bit integer */
typedef char Char;
typedef unsigned char Uchar;
/* default floating type */
typedef double Doub;
typedef long double Ldoub;
typedef bool Bool;
/* semantically ACGT */
#define GF4char Uchar
/* semantically a compressed GF4word */
#define GF4reg Ullong
/* semantically VARIABLE NUMBER of plaintext bits */
#define Mbit Uchar

typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;

#ifdef DPU_ENV
/* host vectors are represented as C Array on DPU side */
typedef Uchar *VecUchar;
/* semantically string of ACGT */
#define GF4word VecUchar
/* message bits unpacked to variable */
#define VecMbit VecUchar
typedef Ullong *VecUllong;
typedef Int *VecInt;
#endif

#ifdef DPU_ENV
#define DPU_HOST_SYMBOL __host
#else
#define DPU_HOST_SYMBOL
#endif

#ifdef DPU_ENV

Int MAX(Int a, Int b) { return b > a ? (b) : (a); }
uint64_t MIN(uint64_t a, uint64_t b) { return b < a ? (b) : (a); }

#endif

#endif
