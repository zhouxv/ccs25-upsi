
#include "./Common.h"

size_t sparse_comp::point_encoding_block_count(size_t d) {
    return (size_t) ceil(((double) d)/4.0);
}