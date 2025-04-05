#!/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ ! -d "${script_dir}/log" ]; then
    mkdir -p "${script_dir}/log"
fi