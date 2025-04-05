#!/bin/bash

echo -e "\e[1;31mAre you sure you want to execute this script? (yes/no): \e[0m\c"
read response

if [[ "$response" != "yes" && "$response" != "y" ]]; then
    echo "Aborting script execution."
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

rm -rf "$SCRIPT_DIR/log"
mkdir -p "$SCRIPT_DIR/log"

rm -rf "$SCRIPT_DIR/bench_results.out"

echo -e "\e[1;32mAll logs and bench-results.out have been deleted.\e[0m"