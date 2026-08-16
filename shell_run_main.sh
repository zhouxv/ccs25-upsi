#!/bin/bash

# Pre-version FPSI benchmark script
#
# LInfPre / L1Pre:
#   d     = 2, 6, 10
#   n     = 8, 12, 16   (set size = 2^n)
#   delta = 10, 60, 250
#
# L2Pre:
#   d     = 2
#   n     = 8, 12, 16
#   delta = 10, 60, 250
#
# metric:
#   m = 0: LInfPre
#   m = 1: L1Pre
#   m = 2: L2Pre

BIN="./build/main"
LOG_FILE="OOTEST_$(date +%Y-%m-%d_%H-%M-%S).log"

# Set default parameters for the benchmark.
trait=20
target_matching_points=29
ms=(0 1 2)
ns=(8 12 16)
dims=(2 6 10)
deltas=(10 60 250)

# ============================================================
# Parse command-line arguments.
# ============================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    -trait)
      trait="$2"
      shift 2
      ;;

    -i)
      target_matching_points="$2"
      shift 2
      ;;

    -m)
      shift
      ms=()
      while [[ $# -gt 0 && "$1" != -* ]]; do
        ms+=("$1")
        shift
      done
      ;;

    -n)
      shift
      ns=()
      while [[ $# -gt 0 && "$1" != -* ]]; do
        ns+=("$1")
        shift
      done
      ;;

    -d)
      shift
      dims=()
      while [[ $# -gt 0 && "$1" != -* ]]; do
        dims+=("$1")
        shift
      done
      ;;

    -delta)
      shift
      deltas=()
      while [[ $# -gt 0 && "$1" != -* ]]; do
        deltas+=("$1")
        shift
      done
      ;;

    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Print the header of the benchmark results.
print_header() {
  printf "%-7s  %7s  %7s  %2s  %7s  %9s  %9s  %9s  %9s  %9s  %9s\n" \
    "" \
    "N" \
    "Metric" \
    "d" \
    "delta" \
    "Off.Com(MB)" \
    "Off.Time(s)" \
    "On.Com(MB)" \
    "On.Time(s)" \
    "Total.Com" \
    "Total.Time" \
    | tee -a "${LOG_FILE}"
}

# ============================================================
# Run benchmarks.
# ============================================================

for metric in "${ms[@]}"; do
  case "${metric}" in
    0)
      echo "==================== LInfPre ====================" | tee -a "${LOG_FILE}"
      print_header

      for dim in "${dims[@]}"; do
        for n in "${ns[@]}"; do
          for delta in "${deltas[@]}"; do
            "${BIN}" \
              -m 0 \
              -n "${n}" \
              -d "${dim}" \
              -delta "${delta}" \
              -i "${target_matching_points}" \
              -trait "${trait}" \
              -log 0 \
              | tee -a "${LOG_FILE}"
          done
        done
      done
      ;;

    1)
      echo "===================== L1Pre =====================" | tee -a "${LOG_FILE}"
      print_header

      for dim in "${dims[@]}"; do
        for n in "${ns[@]}"; do
          for delta in "${deltas[@]}"; do
            "${BIN}" \
              -m 1 \
              -n "${n}" \
              -d "${dim}" \
              -delta "${delta}" \
              -i "${target_matching_points}" \
              -trait "${trait}" \
              -log 0 \
              | tee -a "${LOG_FILE}"
          done
        done
      done
      ;;

    2)
      echo "===================== L2Pre =====================" | tee -a "${LOG_FILE}"
      print_header

      dim=2

      for n in "${ns[@]}"; do
        for delta in "${deltas[@]}"; do
          "${BIN}" \
            -m 2 \
            -n "${n}" \
            -d "${dim}" \
            -delta "${delta}" \
            -i "${target_matching_points}" \
            -trait "${trait}" \
            -log 0 \
            | tee -a "${LOG_FILE}"
        done
      done
      ;;

    *)
      echo "Invalid metric: ${metric}. Supported values: 0, 1, 2."
      exit 1
      ;;
  esac
done

echo "=================================================" | tee -a "${LOG_FILE}"
echo "All benchmarks finished. Results: ${LOG_FILE}" | tee -a "${LOG_FILE}"