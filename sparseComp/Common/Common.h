#pragma once

#include <stdint.h>
#include <algorithm>
#include <array>
#include <type_traits>
#include <cmath>
#include "cryptoTools/Common/block.h"

#define MAX_DIM_DEFINE 10

using osuCrypto::block;

namespace sparse_comp {

    struct point {
        
        static constexpr uint8_t MAX_DIM = MAX_DIM_DEFINE;
        
        uint8_t coord_dim;
        uint32_t coords[MAX_DIM_DEFINE] = {};

        point() {
            this->coord_dim = 0;
        }

        point(uint8_t coord_dim, uint32_t coords[MAX_DIM_DEFINE]) {
            this->coord_dim = (uint8_t) std::min(sparse_comp::point::MAX_DIM, coord_dim);

            for(uint8_t i=0;i < this->coord_dim;i++) {
                this->coords[i] = coords[i];
            }
        }

        point(uint32_t v0, uint32_t v1) {
            this->coord_dim = 2;
            this->coords[0] = v0;
            this->coords[1] = v1;
        }

        uint32_t& operator[](size_t index) {
            return coords[index];
        }

        const uint32_t& operator[](size_t index) const {
            return coords[index];
        }
    };

    template<typename T, size_t n>
    inline void free_array(std::array<T*,n>* arr) {
        for (size_t i=0;i < n;i++) {
            delete arr->at(i);
        }

        delete arr;
    };

    size_t point_encoding_block_count(size_t d);

    template<size_t t, size_t d>
    void points_to_blocks(std::array<point,t>& points, std::vector<block>& blocks) {
        static_assert(d <= point::MAX_DIM);
        const size_t pt_blk_cnt = point_encoding_block_count(d);

        blocks.resize(t*pt_blk_cnt);

        for (size_t i = 0; i < t; i++) {
            uint32_t* coords_ptr = points[i].coords;
            
            for (size_t j=0;j < pt_blk_cnt;j++) {
                uint32_t* data_ptr = (uint32_t*) blocks[i*pt_blk_cnt+j].data();

                memcpy(data_ptr, coords_ptr, sizeof(uint32_t)*4);
                
                coords_ptr += 4;
            }
        }

    }

    
    template<size_t d>
    point blocks_to_point(std::vector<block>& blocks) {
        const size_t pt_blk_cnt = point_encoding_block_count(d);

        point pt;
        pt.coord_dim = d;
        uint32_t* coords_ptr = pt.coords;

        for (size_t j=0;j < pt_blk_cnt;j++) {
            uint32_t* data_ptr = (uint32_t*) blocks[j].data();

            memcpy(coords_ptr, data_ptr, sizeof(uint32_t)*4);
            
            coords_ptr += 4;
        }

        return pt;
    }

    template<size_t t, size_t d>
    void blocks_to_points(std::vector<block>& blocks, std::array<point,t>& points) {
        static_assert(d <= point::MAX_DIM);
        const size_t pt_blk_cnt = point_encoding_block_count(d);
        assert(blocks.size() == t * pt_blk_cnt);

        for (size_t i = 0; i < t; i++) {
            uint32_t* coords_ptr = points[i].coords;
            
            for (size_t j=0;j < pt_blk_cnt;j++) {
                uint32_t* data_ptr = (uint32_t*) blocks[i*pt_blk_cnt+j].data();

                memcpy(coords_ptr, data_ptr, sizeof(uint32_t)*4);
                
                coords_ptr += 4;
            }
        }

    }

};

