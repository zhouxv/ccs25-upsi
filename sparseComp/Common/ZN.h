#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/block.h"
#include "./VecMatrix.h"

using PRNG = osuCrypto::PRNG;
using block = osuCrypto::block;

template <typename T>
using vec = std::vector<T>;

template <typename T, size_t N>
using array = std::array<T, N>;

namespace sparse_comp {

    template <uint64_t N>
    struct alignas(8) ZN {

        static_assert(N <= 1152921504606846976, "N cannot be larger than 2^{60}.");

        std::uint64_t val;

        ZN() {
            this->val = 0;
        }

        ZN(uint64_t new_val) {
            this->val = new_val % N;
        }

        // TODO #1: Move this method to an implementation file.
        // TODO #2: Fix the way this method samples a random element from PRNG! (It is insecure.)
        inline static ZN<N> sample(PRNG& prng) {
            return ZN<N>(prng.get<uint64_t>());
        }

        static VecMatrix<ZN<N>>* sample(PRNG& prng, size_t row_count, size_t column_count) {
            VecMatrix<ZN<N>>* mtx = new VecMatrix<ZN<N>>(row_count, column_count);

            for (size_t i=0;i < row_count;i++) {
                vector<ZN<N>>& row = (*mtx)[i];
                for (size_t j=0;j < column_count;j++) {
                    row[j] = ZN<N>::sample(prng);
                }
            }

            return mtx;
        }

        template<size_t t, size_t k>
        static void sample(PRNG& prng, array<array<ZN<N>,k>,t>& output) {
            
            for (size_t i=0;i < t;i++) {
                array<ZN<N>,k>& row = output[i];
                for (size_t j=0;j < k;j++) {
                    row[j] = ZN<N>(prng.get<uint64_t>());
                }

            }

        }   

        static VecMatrix<ZN<N>>* subVec(VecMatrix<ZN<N>>& mtx, vector<ZN<N>>& v) {
            size_t row_count = mtx.row_count();
            size_t column_count = mtx.col_count();

            VecMatrix<ZN<N>>* result = new VecMatrix<ZN<N>>(row_count,column_count);

            for(size_t i=0;i < row_count;i++) {
                vec<ZN<N>>& mtx_row = mtx[i];
                vec<ZN<N>>& result_row = (*result)[i];
                for(size_t j=0;j < column_count;j++) {
                    result_row[j] = mtx_row[j] - v[i];
                }
            }

            return result;
        }

        template<size_t n>
        static VecMatrix<ZN<N>>* subVec(VecMatrix<ZN<N>>& mtx, array<ZN<N>,n>& v) {
            size_t row_count = mtx.row_count();
            size_t column_count = mtx.col_count();

            VecMatrix<ZN<N>>* result = new VecMatrix<ZN<N>>(row_count,column_count);

            for(size_t i=0;i < row_count;i++) {
                vec<ZN<N>>& mtx_row = mtx[i];
                vec<ZN<N>>& result_row = (*result)[i];
                for(size_t j=0;j < column_count;j++) {
                    result_row[j] = mtx_row[j] - v[i];
                }
            }

            return result;
        }

        template<size_t k, size_t n>
        static VecMatrix<ZN<N>>* subVec(array<array<ZN<N>,n>,k>& mtx, array<ZN<N>,k>& v) {
            VecMatrix<ZN<N>>* result = new VecMatrix<ZN<N>>(k,n);

            for(size_t i=0;i < k;i++) {
                array<ZN<N>,n>& mtx_row = mtx[i];
                vec<ZN<N>>& result_row = (*result)[i];
                for(size_t j=0;j < n;j++) {
                    result_row[j] = mtx_row[j] - v[i];
                }
            }

            return result;
        }

        // Let v be a vector and s be a scalar. This method returns a vector v' where v'[i] = v[i] + s. 
        inline static void vec_scalar_add(const vec<ZN<N>>& v, ZN<N> scalar, vec<ZN<N>>& result_vec) {
            for(uint32_t i=0;i < v.size();i++) {
                result_vec[i] = v[i] + scalar;
            }
        }

        // Let v be a vector and s be a scalar. This method returns a vector v' where v'[i] = v[i] - s. 
        inline static void vec_scalar_sub(const vec<ZN<N>>& v, ZN<N> scalar, vec<ZN<N>>& result_vec) {
            for(uint32_t i=0;i < v.size();i++) {
                result_vec[i] = v[i] - scalar;
            }
        }

        template<size_t l>
        inline static ZN<N> add_vec_components(array<ZN<N>,l>& arr) {
            ZN<N> result;

            for (size_t i=0;i < l;i++) {
                result = result + arr[i];
            }

            return result;
        }

        // bv[i] = (block) zv[i]
        inline static void ZNVec_to_BlockVec(const vec<ZN<N>>& zv, vec<block>& bv) {
            for (uint32_t i=0;i < zv.size();i++) {
                bv[i] = (block) zv;
            }
        }

        template<uint8_t l, size_t n> 
        inline void bit_decomp(array<ZN<2>,n>& v, size_t start_offset) const {
            
            static_assert(l <= 64, "l cannot be larger than 64.");
                        
            for (uint8_t i=0;i < l;i++) v[i + start_offset] = ZN<2>(((this->val) >> i) & 1);

        }

        inline ZN<N> add_inv() {
            ZN<N> result;
            result.val = (N - this->val) % N;
            return result;
        }

        ZN<N> operator+(const ZN<N> other) const {
            ZN<N> result;
            result.val = (this->val + other.val) % N;
            return result;
        }

        ZN<N> operator-(const ZN<N> other) const {
            ZN<N> result;
            result.val = (this->val + (N - other.val)) % N;
            return result;
        }

        uint64_t to_uint64_t() {
            return this->val;
        }

        inline int64_t to_int64_t() {
            return (int64_t)this->val;
        }

        size_t to_size_t() {
            return this->val;
        }

        block to_block() {
            return block(0,this->val);
        }

    };

}