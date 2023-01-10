/*
	Copyright (c) 2022 Marcel Pi Nacy
	
	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:
	
	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.
	
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.
*/

/*
    Martin Burtscher and Paruj Ratanaworabhan's paper:
    https://www.researchgate.net/publication/224323445_FPC_A_High-Speed_Compressor_for_Double-Precision_Floating-Point_Data
*/

#ifndef FPC_INCLUDED
#define FPC_INCLUDED
#include <stdint.h>
#include <stddef.h>

#ifndef FPC_CALL
#define FPC_CALL
#endif

#ifndef FPC_ATTR
#define FPC_ATTR
#endif

#ifdef __cplusplus
#define FPC_EXTERN_C_BEGIN extern "C" {
#define FPC_EXTERN_C_END }
#else
#define FPC_EXTERN_C_BEGIN
#define FPC_EXTERN_C_END
#endif

#define FPC_TABLE_SIZE_DEFAULT 32768

#define FPC_LEAST_FREQUENT_LZBC 4
#define FPC_UPPER_BOUND_METADATA(COUNT) ((size_t)((COUNT) + 1) / 2)
#define FPC_UPPER_BOUND_DATA(COUNT) ((size_t)(COUNT) * 8)
#define FPC_UPPER_BOUND(COUNT) FPC_UPPER_BOUND_METADATA((COUNT)) + FPC_UPPER_BOUND_DATA((COUNT))
#define FPC_DEFAULT_HASH_ARGS { 6, 48, 2, 40 }

#define FPC32_UPPER_BOUND_METADATA(COUNT) ((size_t)((COUNT) + 1) / 2)
#define FPC32_UPPER_BOUND_DATA(COUNT) ((size_t)(COUNT) * 4)
#define FPC32_UPPER_BOUND(COUNT) FPC32_UPPER_BOUND_METADATA((COUNT)) + FPC32_UPPER_BOUND_DATA((COUNT))
#define FPC32_DEFAULT_HASH_ARGS { 1, 22, 4, 23 }

FPC_EXTERN_C_BEGIN
typedef struct fpc_hash_args_t
{
    uint8_t fcm_lshift;
    uint8_t fcm_rshift;
    uint8_t dfcm_lshift;
    uint8_t dfcm_rshift;
} fpc_hash_args_t;

typedef struct fpc_context_t
{
    uint64_t* fcm;
    uint64_t* dfcm;
    // The size, in elements, of the array pointed to by "fcm".
    size_t fcm_size;
    // The size, in elements, of the array pointed to by "dfcm".
    size_t dfcm_size;
    // Seed value.
    double seed;
    // Custom options for the FCM and DFCM hash functions.
    fpc_hash_args_t hash_args;
} fpc_context_t;

typedef struct fpc32_context_t
{
    uint32_t* fcm;
    uint32_t* dfcm;
    // The size, in elements, of the array pointed to by "fcm".
    size_t fcm_size;
    // The size, in elements, of the array pointed to by "dfcm".
    size_t dfcm_size;
    // Seed value.
    float seed;
    // Custom options for the FCM and DFCM hash functions.
    fpc_hash_args_t hash_args;
} fpc32_context_t;

FPC_ATTR size_t FPC_CALL fpc_encode_explicit(
    fpc_context_t* ctx,
    const double* values,
    size_t value_count,
    void* out_headers,
    void* out_compressed);

FPC_ATTR void FPC_CALL fpc_decode_explicit(
    fpc_context_t* ctx,
    const void* headers,
    const void* compressed,
    double* out_values,
    size_t out_count);

FPC_ATTR size_t FPC_CALL fpc_encode(
    fpc_context_t* ctx,
    const double* values,
    size_t value_count,
    void* out_compressed);

FPC_ATTR void FPC_CALL fpc_decode(
    fpc_context_t* ctx,
    const void* compressed,
    double* out_values,
    size_t out_count);

FPC_ATTR size_t FPC_CALL fpc32_encode_explicit(
    fpc32_context_t* ctx,
    const float* values,
    size_t value_count,
    void* out_headers,
    void* out_compressed);

FPC_ATTR void FPC_CALL fpc32_decode_explicit(
    fpc32_context_t* ctx,
    const void* headers,
    const void* compressed,
    float* out_values,
    size_t out_count);

FPC_ATTR size_t FPC_CALL fpc32_encode(
    fpc32_context_t* ctx,
    const float* values,
    size_t value_count,
    void* out_compressed);

FPC_ATTR void FPC_CALL fpc32_decode(
    fpc32_context_t* ctx,
    const void* compressed,
    float* out_values,
    size_t out_count);
FPC_EXTERN_C_END
#endif