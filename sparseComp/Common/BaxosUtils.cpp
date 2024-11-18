#include "./BaxosUtils.h"
#include <cmath>

uint64_t sparse_comp::baxosBinSize(size_t itemCount) {
    return (uint64_t) ceil(1.27*((double) itemCount));
}