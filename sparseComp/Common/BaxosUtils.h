#pragma once
#include <cstdint>
#include <cstddef>

namespace sparse_comp {

    uint64_t baxosBinSize(size_t itemCount);
    size_t baxosBlockCount(size_t itemCount, size_t ssp);

};