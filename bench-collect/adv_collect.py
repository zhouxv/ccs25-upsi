import subprocess
import re

OKVS_BENCH_EXE_PATH = "../build/okvs_bench"
ENC_BENCH_EXE_PATH = "../build/enc_bench"
OPRF_BENCH_EXE_PATH = "../build/oprf_bench"

BENCH_PARAMS = [
    {"log2n": 8, "d": 2, "delta": 10},
    #{"log2n": 8, "d": 6, "delta": 10},
    #{"log2n": 8, "d": 10, "delta": 10},
    #{"log2n": 8, "d": 2, "delta": 30},
    #{"log2n": 8, "d": 6, "delta": 30},
    #{"log2n": 8, "d": 10, "delta": 30},
    #{"log2n": 12, "d": 2, "delta": 10},
    #{"log2n": 12, "d": 6, "delta": 10},
    #{"log2n": 12, "d": 10, "delta": 10},
    #{"log2n": 12, "d": 2, "delta": 30},
    #{"log2n": 12, "d": 6, "delta": 30},
    #{"log2n": 12, "d": 10, "delta": 30},
    #{"log2n": 16, "d": 2, "delta": 10},
    #{"log2n": 16, "d": 6, "delta": 10},
    #{"log2n": 16, "d": 10, "delta": 10},
    #{"log2n": 16, "d": 2, "delta": 30},
    #{"log2n": 16, "d": 6, "delta": 30},
    #{"log2n": 16, "d": 10, "delta": 30}
]

def extract_oprf_measurement(cmd_str_out):
    match = re.search(r'(\d+(\.\d+)?)\s*ms', cmd_str_out)
    time_ms = float(match.group(1)) if match else None
    
    match_str = re.search(r'Number of MBs exchanged:\s+(\d+(\.\d+)?)', cmd_str_out)
    number_after_str = float(match_str.group(1)) if match_str else None
    return time_ms, number_after_str

def extract_aes_measurement(cmd_str_out):
    match = re.search(r'(\d+(\.\d+)?)\s*ms', cmd_str_out)
    return float(match.group(1)) if match else None
    
for params in BENCH_PARAMS:
    log2n = params["log2n"]
    d = params["d"]
    delta = params["delta"]
    
    # cmd = [OPRF_BENCH_EXE_PATH, "--benchmark-samples", "1", f"[oprf][nA=nB=2^{log2n}][d={d}]"]
    cmd = [ENC_BENCH_EXE_PATH, "--benchmark-samples", "1", f"[enc_dec][nA=nB=2^{log2n}][d={d}]"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    o, e = proc.communicate()

    print(o.decode())
    print(extract_aes_measurement(o.decode()))

    # print(o.decode())
    # print(extract_oprf_measurement(o.decode()))