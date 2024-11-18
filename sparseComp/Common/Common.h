#pragma once

#include <stdint.h>
#include <algorithm>
#include <array>
#include <type_traits>


#define MAX_DIM_DEFINE 10

namespace sparse_comp {

    struct point {
        
        static const uint8_t MAX_DIM;
        
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
    };

    template<typename T, size_t n>
    inline void free_array(std::array<T*,n>* arr) {
        for (size_t i=0;i < n;i++) {
            delete arr->at(i);
        }

        delete arr;
    }

};

