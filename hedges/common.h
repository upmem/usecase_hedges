#ifndef __COMMON__
#define __COMMON__
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#ifndef DPU_ENV
#include <vector>
#else
#include <mram.h>
#include <alloc.h>
#endif

/* Size of the buffer for which we compute the checksum: 64KBytes. */
#define BUFFER_SIZE (1 << BUFFER_SIZE_LOG2)

// basic type names (redefine if your bit lengths don't match)

typedef int Int; // 32 bit integer
typedef unsigned int Uint;
typedef char Char; // 8 bit integer
typedef unsigned char Uchar;
typedef double Doub; // default floating type
typedef long double Ldoub;
typedef bool Bool;

#define GF4char Uchar // semantically ACGT
#define GF4reg Ullong // semantically a compressed GF4word
#define Mbit Uchar    // semantically VARIABLE NUMBER of plaintext bits

typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;

#ifdef DPU_ENV
// host vectors are represented as C Array on DPU side
typedef Uchar *VecUchar;
#define GF4word VecUchar // semantically string of ACGT
#define VecMbit VecUchar // message bits unpacked to variable
typedef Ullong *VecUllong;
typedef Int *VecInt;
#endif

#ifdef DPU_ENV
#define DPU_HOST_SYMBOL __host
#else
#define DPU_HOST_SYMBOL
#endif

/**
 *  extra Vector Size DPU/HOST Shared varibles
 *  required for SHARED DPU/HOST vector codding stype (array dynamicly allocated on both DPU/HOST)
 **/

#ifdef DPU_ENV

/**
 *  Scalar DPU/HOST variables
 * **/

/**
 template <class T>
 inline const T &MAX(const T &a, const T &b)
 {
     return b > a ? (b) : (a);
}

float MAX(double a, float b)
{
    return b > a ? (b) : float(a);
}

float MAX(float a, double b)
{
    return b > a ? float(b) : (a);
}
**/
Int MAX(Int a, Int b)
{
    return b > a ? (b) : (a);
}
uint64_t MIN(uint64_t a, uint64_t b)
{
    return b < a ? (b) : (a);
}

#endif

#endif /* -_COMMON_- */
