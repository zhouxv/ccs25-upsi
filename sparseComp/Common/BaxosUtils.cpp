#include "./BaxosUtils.h"
#include "volePSI/Paxos.h"
#include <cmath>

using Baxos = volePSI::Baxos;

uint64_t sparse_comp::baxosBinSize(size_t itemCount) {
    return (uint64_t) ceil(1.27*((double) itemCount));
}

size_t sparse_comp::baxosBlockCount(size_t itemCount, size_t ssp) {

    Baxos paxos;
    paxos.init(itemCount, sparse_comp::baxosBinSize(itemCount), 3, ssp, volePSI::PaxosParam::GF128, oc::ZeroBlock);

    return paxos.size();

}