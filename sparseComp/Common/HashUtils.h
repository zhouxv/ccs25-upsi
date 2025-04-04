#pragma once

#include "./Common.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/AES.h"
#include <array>

using block = osuCrypto::block;
using AES = osuCrypto::AES;
using point = sparse_comp::point;

static constexpr size_t int_pow(size_t base, size_t exp) {
    return (exp == 0) ? 1 : base * int_pow(base, exp - 1);
}

namespace sparse_comp {

    block hash_point(const AES& aes, const point& point);
    block hash_point(const AES& aes, const point& point, size_t sot_idx, size_t msg_vec_idx);
    
    inline point spatial_hash(const point& in_point, size_t point_dim, uint8_t delta) {
        point out_point;
        out_point.coord_dim = point_dim;

        for (size_t j = 0; j < point_dim; j++) {
            out_point.coords[j] = in_point.coords[j] / ((uint32_t) (2*delta));
        }

        return out_point;
    }

    template<size_t t>
    inline void spatial_hash(std::array<point,t>& in_points, std::array<point,t>& out_points, size_t point_dim, uint8_t delta) {
        for (size_t i = 0; i < t; i++) {
            out_points[i].coord_dim = point_dim;

            for (size_t j = 0; j < point_dim; j++) {
                out_points[i].coords[j] = in_points[i].coords[j] / ((uint32_t) (2*delta));
            }
        }
    }

    template<size_t t, size_t d, size_t cell_count>
    inline void spatial_cell_hash(std::array<point, t>& point_center, std::array<point, cell_count>& cells, uint8_t delta) {
        constexpr const size_t twotod = int_pow(2, d);

        static_assert(cell_count == twotod*t);
        static_assert(d <= point::MAX_DIM);

        uint32_t u32_delta = ((uint32_t) delta); 
        
        for (size_t i = 0; i < t; i++) {
            for (size_t j = 0; j < twotod; j++) {
                
                cells[twotod*i+j].coord_dim = point_center[i].coord_dim;

                for (size_t k = 0; k < d; k++) {
                    
                    cells[twotod*i+j].coords[k] = point_center[i].coords[k]/(2*u32_delta);
                    
                    uint32_t b = ((j >> k) & 1);
                    if (b == 0) continue;

                    if((point_center[i].coords[k] + u32_delta)/(2*u32_delta) > (point_center[i].coords[k])/(2*u32_delta)) {
                        cells[twotod*i+j].coords[k] += b;
                    } else if(b == 1) {
                        cells[twotod*i+j].coords[k] -= b;
                    }
                }

            }
        }

    }

};