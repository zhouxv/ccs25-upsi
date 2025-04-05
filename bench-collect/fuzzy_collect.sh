#!/bin/bash

N_BENCH_SAMPLES=1
METRIC="fuzzylinf"
TEST_BIN_NAME="${METRIC}_bench"
OUT_FILE_PATH="bench_results.out"

BENCH_PARAMS=(
    "256 2 10 40"
    "256 6 10 40"
    #"256 10 10 40"
    #"256 2 30 40"
    #"256 6 30 40"
    #"256 10 30 40"
    #"4096 2 10 40"
    #"4096 6 10 40"
    #"4096 10 10 40"
    #"4096 2 30 40"
    #"4096 6 30 40"
    #"4096 10 30 40"
    #"65536 2 10 40"
    #"65536 6 10 40"
    #"65536 10 10 40"
    #"65536 2 30 40"
    #"65536 6 30 40"
    #"65536 10 30 40"
)

SCRIPT_PATH=$(realpath "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")

BUILD_PATH="${SCRIPT_DIR}/../build"
LOG_PATH="${SCRIPT_DIR}/log"
TEST_BIN_PATH="${BUILD_PATH}/${TEST_BIN_NAME}"

test_name() {
    N=$1
    D=$2
    DELTA=$3
    SSP=$4
    METRIC=$5

    echo "${METRIC} (n=m=${N} d=${D} delta=${DELTA} ssp=${SSP})"
}

extract_metrics() {
    BENCH_OUTPUT=$1

    LAST_BENCH_NAME_LINDEX=$(echo "$BENCH_OUTPUT" | sed -n "/benchmark name/=" | tail -n 1)
    RUNTIME_LINE=$(echo "$BENCH_OUTPUT" | tail -n +"$LAST_BENCH_NAME_LINDEX" | head -n 6 | tail -n 2 | head -n 1)

    RUNTIME=$(echo "$RUNTIME_LINE" | awk '{print $1, $2}')

    #echo "RUNTIME: $RUNTIME"
    MBs_EXCHANGED_LINE=$(echo "$BENCH_OUTPUT" | grep "Number of MBs exchanged")
    MBs_EXCHANGED=$(echo "$MBs_EXCHANGED_LINE" | grep -oP '\d+\.\d+' | tail -n 1)

    echo "$RUNTIME;$MBs_EXCHANGED"

}

RUNTIME_OUT_STR=""
MBs_EXCHANGED_OUT_STR=""

for PARAMS in "${BENCH_PARAMS[@]}"; do
    set -- $PARAMS
    
    TEST_NAME=$(test_name "$1" "$2" "$3" "$4" "$METRIC")
    
    echo "Running test: $TEST_NAME"
    
    BENCH_OUTPUT=$("$TEST_BIN_PATH" --benchmark-samples $N_BENCH_SAMPLES --success "$TEST_NAME")
    #echo $BENCH_OUTPUT
    echo $BENCH_OUTPUT > "${LOG_PATH}/${TEST_NAME}.out"
    METRICS=$(extract_metrics "$BENCH_OUTPUT")
    
    MBs_EXCHANGED=$(echo "$METRICS" | cut -d ';' -f 2)
    RUNTIME=$(echo "$METRICS" | cut -d ';' -f 1)

    echo ""
    echo "RUNTIME: $RUNTIME"
    echo "MBs exchanged: $MBs_EXCHANGED"
    echo ""

    RUNTIME_OUT_STR="${RUNTIME_OUT_STR}${RUNTIME}\t"
    MBs_EXCHANGED_OUT_STR="${MBs_EXCHANGED_OUT_STR}${MBs_EXCHANGED}\t"
done

printf "$RUNTIME_OUT_STR\n" > "$OUT_FILE_PATH"
printf "$MBs_EXCHANGED_OUT_STR\n" >> "$OUT_FILE_PATH"
