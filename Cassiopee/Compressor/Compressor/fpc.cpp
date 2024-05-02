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

#include "fpc.h"

//================================================================
// FPC Implementation:
//================================================================

#if defined(__clang__) || defined(__GNUC__)
  #define FPC_GCC_OR_CLANG
  #define FPC_ALIGN(K) __attribute__((aligned((K))))
  #define FPC_CLZ32(MASK) (uint_fast8_t)__builtin_clz(MASK)
  #define FPC_CLZ64(MASK) (uint_fast8_t)__builtin_clzll(MASK)

#elif defined(_MSC_VER) || defined(_MSVC_LANG)
  #include <intrin.h>
  #include <Windows.h>
  #define FPC_MSVC
  #define FPC_ALIGN(K) __declspec(align(K))

  #if defined(_M_IX86) || defined(_M_X64)
    #define FPC_CLZ32(MASK) (uint_fast8_t)__lzcnt(MASK)
    #define FPC_CLZ64(MASK) (uint_fast8_t)__lzcnt64(MASK)
  #elif defined(_M_ARM)
    #define FPC_CLZ32(MASK) (uint_fast8_t)_CountLeadingZeros(MASK)
    #define FPC_CLZ64(MASK) (uint_fast8_t)_CountLeadingZeros64(MASK)
  #else
    #error "FPC: UNSUPPORTED ARCHITECTURE."
  #endif

#else
  #error "FPC: UNSUPPORTED COMPILER"
#endif

#if defined(_DEBUG) || !defined(NDEBUG)
  #include <assert.h>
  #define FPC_INVARIANT(EXPRESSION) assert((EXPRESSION))
#else
  #define FPC_INVARIANT(EXPRESSION) 
  //#ifdef FPC_GCC_OR_CLANG
  //  #define FPC_INVARIANT(EXPRESSION) __builtin_assume((EXPRESSION))
  //#elif defined(FPC_MSVC)
  //  #define FPC_INVARIANT(EXPRESSION) __assume((EXPRESSION))
  //#endif
#endif

#define FPC_FCM_HASH_UPDATE(HASH, VALUE)   HASH = ((HASH << ctx->hash_args.fcm_lshift) ^ (size_t)((VALUE) >> ctx->hash_args.fcm_rshift)) & fcm_mod_mask
#define FPC_DFCM_HASH_UPDATE(HASH, VALUE)  HASH = ((HASH << ctx->hash_args.dfcm_lshift) ^ (size_t)((VALUE) >> ctx->hash_args.dfcm_rshift)) & dfcm_mod_mask

#define FPC_IS_ALIGNED(PTR, ALIGN) (((size_t)(PTR) & (size_t)((ALIGN) - 1)) == 0)

#include <string.h>

FPC_EXTERN_C_BEGIN
FPC_ATTR size_t FPC_CALL fpc_encode_explicit(
    fpc_context_t* ctx,
    const double* values,
    size_t value_count,
    void* out_headers,
    void* out_compressed)
{
    const double* end;
    uint8_t* compressed_begin;
    uint8_t* out;
    uint8_t* out_h;
    uint8_t type_a, type_b, lzbc_a, lzbc_b;
    size_t fcm_mod_mask, dfcm_mod_mask, fcm_hash, dfcm_hash;
    uint64_t fcm_prediction, dfcm_prediction, dfcm_delta, prior_value, fcm_xor, dfcm_xor;
    FPC_ALIGN(16)
    uint64_t tmp[2];
    FPC_INVARIANT((uint8_t*)out_headers != (uint8_t*)out_compressed);
    FPC_INVARIANT(
        (uint8_t*)out_headers < (uint8_t*)out_compressed ||
        (uint8_t*)out_headers + value_count >= (uint8_t*)out_compressed);
    out_h = (uint8_t*)out_headers;
    end = values + value_count;
    compressed_begin = out = (uint8_t*)out_compressed;
    prior_value = *(const uint64_t*)&ctx->seed;
    fcm_prediction = dfcm_prediction = fcm_hash = dfcm_hash = 0;
    fcm_mod_mask = ctx->fcm_size - 1;
    dfcm_mod_mask = ctx->dfcm_size - 1;
    while ((end - values) >= 2)
    {
        (void)memcpy(tmp, values, 16);
        values += 2;
        // FCM STEP #1
        fcm_xor = tmp[0] ^ fcm_prediction;
        ctx->fcm[fcm_hash] = tmp[0];
        FPC_FCM_HASH_UPDATE(fcm_hash, tmp[0]);
        fcm_prediction = ctx->fcm[fcm_hash];
        // DFCM STEP #1
        dfcm_delta = tmp[0] - prior_value;
        dfcm_xor = tmp[0] ^ (dfcm_prediction + prior_value);
        prior_value = tmp[0];
        ctx->dfcm[dfcm_hash] = dfcm_delta;
        FPC_DFCM_HASH_UPDATE(dfcm_hash, dfcm_delta);
        dfcm_prediction = ctx->dfcm[dfcm_hash];
        type_a = fcm_xor > dfcm_xor;
        tmp[0] = type_a ? dfcm_xor : fcm_xor;
        type_a <<= 3;
        // FCM STEP #2
        fcm_xor = tmp[1] ^ fcm_prediction;
        ctx->fcm[fcm_hash] = tmp[1];
        FPC_FCM_HASH_UPDATE(fcm_hash, tmp[1]);
        fcm_prediction = ctx->fcm[fcm_hash];
        // DFCM STEP #2
        dfcm_delta = tmp[1] - prior_value;
        dfcm_xor = tmp[1] ^ (dfcm_prediction + prior_value);
        prior_value = tmp[1];
        ctx->dfcm[dfcm_hash] = dfcm_delta;
        FPC_DFCM_HASH_UPDATE(dfcm_hash, dfcm_delta);
        dfcm_prediction = ctx->dfcm[dfcm_hash];
        type_b = fcm_xor > dfcm_xor;
        tmp[1] = type_b ? dfcm_xor : fcm_xor;
        type_b <<= 3;
        // COMPRESS PAIR OF VALUES
        lzbc_a = FPC_CLZ64(tmp[0]) >> 3;
        lzbc_b = FPC_CLZ64(tmp[1]) >> 3;
        lzbc_a -= (lzbc_a == FPC_LEAST_FREQUENT_LZBC);
        lzbc_b -= (lzbc_b == FPC_LEAST_FREQUENT_LZBC);
        *out_h =
            ((type_b | (lzbc_b - (lzbc_b > FPC_LEAST_FREQUENT_LZBC))) << 4) |
            (type_a | (lzbc_a - (lzbc_a > FPC_LEAST_FREQUENT_LZBC)));
        ++out_h;
        lzbc_a = 8 - lzbc_a;
        lzbc_b = 8 - lzbc_b;
        (void)memcpy(out, tmp, lzbc_a);
        out += lzbc_a;
        (void)memcpy(out, tmp + 1, lzbc_b);
        out += lzbc_b;
    }
    if (values != end)
    {
        (void)memcpy(tmp, values, 8);
        // LAST FCM STEP
        fcm_xor = tmp[0] ^ fcm_prediction;
        // LAST DFCM STEP
        dfcm_delta = tmp[0] - prior_value;
        dfcm_xor = tmp[0] ^ (prior_value + dfcm_prediction);
        type_a = fcm_xor > dfcm_xor;
        tmp[0] = type_a ? dfcm_xor : fcm_xor;
        // COMPRESS LAST VALUE
        lzbc_a = FPC_CLZ64(tmp[0]) >> 3;
        lzbc_a -= (lzbc_a == FPC_LEAST_FREQUENT_LZBC);
        *out_h = (type_a << 3) | (lzbc_a - (lzbc_a > FPC_LEAST_FREQUENT_LZBC));
        lzbc_a = 8 - lzbc_a;
        (void)memcpy(out, tmp, lzbc_a);
        out += lzbc_a;
    }
    return (out - compressed_begin) + (value_count + 1) / 2;
}

FPC_ATTR void FPC_CALL fpc_decode_explicit(
    fpc_context_t* ctx,
    const void* headers,
    const void* compressed,
    double* out_values,
    size_t out_count)
{
    uint8_t* in;
    uint8_t* in_headers;
    double* end;
    uint8_t lzbc, header;
    size_t fcm_mod_mask, dfcm_mod_mask, fcm_hash, dfcm_hash;
    uint64_t fcm_prediction, dfcm_prediction, dfcm_delta, prior_value;
    FPC_ALIGN(16)
    uint64_t tmp[2];
    in = (uint8_t*)compressed;
    in_headers = (uint8_t*)headers;
    end = out_values + out_count;
    fcm_mod_mask = ctx->fcm_size - 1;
    dfcm_mod_mask = ctx->dfcm_size - 1;
    (void)memcpy(&prior_value, &ctx->seed, 8);
    fcm_hash = dfcm_hash = 0;
    fcm_prediction = dfcm_prediction = 0;
    while (out_values != end)
    {
        // READ HEADER PAIR
        header = *in_headers;
        ++in_headers;
        // DECODE HEADER #1
        lzbc = (header & 7);
        lzbc += (lzbc >= FPC_LEAST_FREQUENT_LZBC);
        lzbc = 8 - lzbc;
        // READ COMPRESSED VALUE #1
        tmp[0] = 0;
        (void)memcpy(tmp, in, lzbc);
        in += lzbc;
        // DECODE VALUE #1
        tmp[0] ^= (header & 8) ? dfcm_prediction : fcm_prediction;
        // FCM UPDATE #1
        ctx->fcm[fcm_hash] = tmp[0];
        FPC_FCM_HASH_UPDATE(fcm_hash, tmp[0]);
        fcm_prediction = ctx->fcm[fcm_hash];
        // DFCM UPDATE #1
        dfcm_delta = tmp[0] - prior_value;
        ctx->dfcm[dfcm_hash] = dfcm_delta;
        FPC_DFCM_HASH_UPDATE(dfcm_hash, dfcm_delta);
        dfcm_prediction = ctx->dfcm[dfcm_hash] + tmp[0];
        prior_value = tmp[0];
        ++out_values;
        if (out_values == end)
        {
            (void)memcpy(out_values-1, tmp, 8);
            return;
        }
        // DECODE HEADER #2
        header >>= 4;
        lzbc = (header & 7);
        lzbc += (lzbc >= FPC_LEAST_FREQUENT_LZBC);
        lzbc = 8 - lzbc;
        // READ COMPRESSED VALUE #2
        tmp[1] = 0;
        (void)memcpy(tmp + 1, in, lzbc);
        in += lzbc;
        // DECODE VALUE #2
        tmp[1] ^= (header & 8) ? dfcm_prediction : fcm_prediction;
        // FCM UPDATE #2
        ctx->fcm[fcm_hash] = tmp[1];
        FPC_FCM_HASH_UPDATE(fcm_hash, tmp[1]);
        fcm_prediction = ctx->fcm[fcm_hash];
        // DFCM UPDATE #2
        dfcm_delta = tmp[1] - prior_value;
        ctx->dfcm[dfcm_hash] = dfcm_delta;
        FPC_DFCM_HASH_UPDATE(dfcm_hash, dfcm_delta);
        dfcm_prediction = ctx->dfcm[dfcm_hash] + tmp[1];
        prior_value = tmp[1];
        // STORE DECOMPRESSED PAIR
        (void)memcpy(out_values - 1, tmp, 16);
        ++out_values;
    }
}

FPC_ATTR size_t FPC_CALL fpc_encode(
    fpc_context_t* ctx,
    const double* values,
    size_t value_count,
    void* out_compressed)
{
    return fpc_encode_explicit(ctx, values, value_count, out_compressed, (uint8_t*)out_compressed + FPC_UPPER_BOUND_METADATA(value_count));
}

FPC_ATTR void FPC_CALL fpc_decode(
    fpc_context_t* ctx,
    const void* compressed,
    double* out_values,
    size_t out_count)
{
    fpc_decode_explicit(ctx, compressed, (uint8_t*)compressed + FPC_UPPER_BOUND_METADATA(out_count), out_values, out_count);
}

FPC_ATTR size_t FPC_CALL fpc32_encode_explicit(
    fpc32_context_t* ctx,
    const float* values,
    size_t value_count,
    void* out_headers,
    void* out_compressed)
{
    const float* end;
    uint8_t* out;
    uint8_t* out_h;
    uint8_t type_a, type_b, lzbc_a, lzbc_b;
    size_t fcm_mod_mask, dfcm_mod_mask, fcm_hash, dfcm_hash;
    uint32_t fcm_prediction, dfcm_prediction, dfcm_delta, prior_value, fcm_xor, dfcm_xor;
    FPC_ALIGN(8)
    uint32_t tmp[2];
    FPC_INVARIANT((uint8_t*)out_headers != (uint8_t*)out_compressed);
    FPC_INVARIANT(
        (uint8_t*)out_headers < (uint8_t*)out_compressed ||
        (uint8_t*)out_headers + value_count >= (uint8_t*)out_compressed);
    end = values + value_count;
    out = (uint8_t*)out_compressed;
    out_h = (uint8_t*)out_headers;
    prior_value = *(const uint32_t*)&ctx->seed;
    fcm_prediction = dfcm_prediction = 0;
    fcm_hash = dfcm_hash = 0;
    fcm_mod_mask = ctx->fcm_size - 1;
    dfcm_mod_mask = ctx->dfcm_size - 1;
    while ((end - values) >= 2)
    {
        (void)memcpy(tmp, values, 8);
        values += 2;
        // FCM STEP #1
        fcm_xor = tmp[0] ^ fcm_prediction;
        ctx->fcm[fcm_hash] = tmp[0];
        FPC_FCM_HASH_UPDATE(fcm_hash, tmp[0]);
        fcm_prediction = ctx->fcm[fcm_hash];
        // DFCM STEP #1
        dfcm_delta = tmp[0] - prior_value;
        dfcm_xor = tmp[0] ^ (dfcm_prediction + prior_value);
        prior_value = tmp[0];
        ctx->dfcm[dfcm_hash] = dfcm_delta;
        FPC_DFCM_HASH_UPDATE(dfcm_hash, dfcm_delta);
        dfcm_prediction = ctx->dfcm[dfcm_hash];
        type_a = fcm_xor > dfcm_xor;
        tmp[0] = type_a ? dfcm_xor : fcm_xor;
        // FCM STEP #2
        fcm_xor = tmp[1] ^ fcm_prediction;
        ctx->fcm[fcm_hash] = tmp[1];
        FPC_FCM_HASH_UPDATE(fcm_hash, tmp[1]);
        fcm_prediction = ctx->fcm[fcm_hash];
        // DFCM STEP #2
        dfcm_delta = tmp[1] - prior_value;
        dfcm_xor = tmp[1] ^ (dfcm_prediction + prior_value);
        prior_value = tmp[1];
        ctx->dfcm[dfcm_hash] = dfcm_delta;
        FPC_DFCM_HASH_UPDATE(dfcm_hash, dfcm_delta);
        dfcm_prediction = ctx->dfcm[dfcm_hash];
        type_b = fcm_xor > dfcm_xor;
        tmp[1] = type_b ? dfcm_xor : fcm_xor;
        // COMPRESS PAIR OF VALUES
        lzbc_a = FPC_CLZ32(tmp[0]) >> 3;
        lzbc_b = FPC_CLZ32(tmp[1]) >> 3;
        type_a <<= 3;
        type_b <<= 3;
        *out_h = ((type_b | lzbc_b) << 4) | (type_a | lzbc_a);
        ++out_h;
        lzbc_a = 4 - lzbc_a;
        lzbc_b = 4 - lzbc_b;
        (void)memcpy(out, tmp, lzbc_a);
        out += lzbc_a;
        (void)memcpy(out, tmp + 1, lzbc_b);
        out += lzbc_b;
    }
    if (values != end)
    {
        (void)memcpy(tmp, values, 4);
        // LAST FCM STEP
        fcm_xor = tmp[0] ^ fcm_prediction;
        // LAST DFCM STEP
        dfcm_delta = tmp[0] - prior_value;
        dfcm_xor = tmp[0] ^ (prior_value + dfcm_prediction);
        type_a = fcm_xor > dfcm_xor;
        tmp[0] = type_a ? dfcm_xor : fcm_xor;
        // COMPRESS LAST VALUE
        lzbc_a = FPC_CLZ32(tmp[0]) >> 3;
        lzbc_a -= (lzbc_a == FPC_LEAST_FREQUENT_LZBC);
        *out_h = (type_a << 3) | (lzbc_a - (lzbc_a > FPC_LEAST_FREQUENT_LZBC));
        lzbc_a = 4 - lzbc_a;
        (void)memcpy(out, tmp, lzbc_a);
        out += lzbc_a;
    }
    return (out - (uint8_t*)out_compressed) + (value_count + 1) / 2;
}

FPC_ATTR void FPC_CALL fpc32_decode_explicit(
    fpc32_context_t* ctx,
    const void* headers,
    const void* compressed,
    float* out_values,
    size_t out_count)
{
    uint8_t* in;
    uint8_t* in_headers;
    float* end;
    uint8_t lzbc, header;
    size_t fcm_mod_mask, dfcm_mod_mask, fcm_hash, dfcm_hash;
    uint32_t fcm_prediction, dfcm_prediction, dfcm_delta, prior_value;
    FPC_ALIGN(8)
    uint32_t tmp[2];
    in = (uint8_t*)compressed;
    in_headers = (uint8_t*)headers;
    end = out_values + out_count;
    fcm_mod_mask = ctx->fcm_size - 1;
    dfcm_mod_mask = ctx->dfcm_size - 1;
    (void)memcpy(&prior_value, &ctx->seed, 4);
    fcm_hash = dfcm_hash = 0;
    fcm_prediction = dfcm_prediction = 0;
    while (out_values != end)
    {
        // READ HEADER PAIR
        header = *in_headers;
        ++in_headers;
        // DECODE HEADER #1
        lzbc = (header & 7);
        lzbc = 4 - lzbc;
        // READ COMPRESSED VALUE #1
        tmp[0] = 0;
        (void)memcpy(tmp, in, lzbc);
        in += lzbc;
        // DECODE VALUE #1
        tmp[0] ^= (header & 8) ? dfcm_prediction : fcm_prediction;
        // FCM UPDATE #1
        ctx->fcm[fcm_hash] = tmp[0];
        FPC_FCM_HASH_UPDATE(fcm_hash, tmp[0]);
        fcm_prediction = ctx->fcm[fcm_hash];
        // DFCM UPDATE #1
        dfcm_delta = tmp[0] - prior_value;
        ctx->dfcm[dfcm_hash] = dfcm_delta;
        FPC_DFCM_HASH_UPDATE(dfcm_hash, dfcm_delta);
        dfcm_prediction = ctx->dfcm[dfcm_hash] + tmp[0];
        prior_value = tmp[0];
        ++out_values;
        if (out_values == end)
        {
            (void)memcpy(out_values-1, tmp, 4);
            return;
        }
        // DECODE HEADER #2
        header >>= 4;
        lzbc = (header & 7);
        lzbc = 4 - lzbc;
        // READ COMPRESSED VALUE #2
        tmp[1] = 0;
        (void)memcpy(tmp + 1, in, lzbc);
        in += lzbc;
        // DECODE VALUE #2
        tmp[1] ^= (header & 8) ? dfcm_prediction : fcm_prediction;
        // FCM UPDATE #2
        ctx->fcm[fcm_hash] = tmp[1];
        FPC_FCM_HASH_UPDATE(fcm_hash, tmp[1]);
        fcm_prediction = ctx->fcm[fcm_hash];
        // DFCM UPDATE #2
        dfcm_delta = tmp[1] - prior_value;
        ctx->dfcm[dfcm_hash] = dfcm_delta;
        FPC_DFCM_HASH_UPDATE(dfcm_hash, dfcm_delta);
        dfcm_prediction = ctx->dfcm[dfcm_hash] + tmp[1];
        prior_value = tmp[1];
        // STORE DECOMPRESSED PAIR
        (void)memcpy(out_values - 1, tmp, 8);
        ++out_values;
    }
}

FPC_ATTR size_t FPC_CALL fpc32_encode(
    fpc32_context_t* ctx,
    const float* values,
    size_t value_count,
    void* out_compressed)
{
    return fpc32_encode_explicit(ctx, values, value_count, out_compressed, (uint8_t*)out_compressed + FPC32_UPPER_BOUND_METADATA(value_count));
}

FPC_ATTR void FPC_CALL fpc32_decode(
    fpc32_context_t* ctx,
    const void* compressed,
    float* out_values,
    size_t out_count)
{
    fpc32_decode_explicit(ctx, compressed, (uint8_t*)compressed + FPC32_UPPER_BOUND_METADATA(out_count), out_values, out_count);
}
FPC_EXTERN_C_END
