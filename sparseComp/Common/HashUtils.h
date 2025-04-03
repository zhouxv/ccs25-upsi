#pragma once

#include "./Common.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/AES.h"
#include <array>

using block = osuCrypto::block;
using AES = osuCrypto::AES;
using point = sparse_comp::point;

namespace sparse_comp {

    block hash_point(const AES& aes, const point& point);
    block hash_point(const AES& aes, const point& point, size_t sot_idx, size_t msg_vec_idx);
    
    template<size_t t>
    inline void spatial_hash(std::array<point,t>& in_points, std::array<point,t>& out_points, size_t point_dim, uint8_t delta) {
        for (size_t i = 0; i < t; i++) {
            out_points[i].coord_dim = point_dim;

            for (size_t j = 0; j < point_dim; j++) {
                out_points[i].coords[j] = in_points[i].coords[j] / (2*delta);
            }
        }
    }

};