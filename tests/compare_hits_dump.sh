#!/usr/bin/env bash
set -uo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
prot="$1" dna="$2"
shift 2
frames=() bin="$repo/build/cuda_alignment"
while [ $# -gt 0 ]; do
    if [ "$1" = -- ]; then bin="$2"; shift 2; else frames+=("$1"); shift; fi
done
[ ${#frames[@]} -gt 0 ] || frames=(3 6)
[ -x "$bin" ] || { echo "No executable at $bin" >&2; exit 2; }

out="$(mktemp -d)"
trap 'rm -rf "$out"' EXIT

fail=0
for frame in "${frames[@]}"; do
    for mode in 0 1; do
        printf '%s\n' "$prot" "$dna" "$mode" "$frame" 0 |
            HITS_DUMP="$out/m$mode.f$frame" "$bin" > "$out/m$mode.f$frame.log" 2>&1 ||
            { echo "mode $mode, $frame frames: exited with $?"; tail -5 "$out/m$mode.f$frame.log"; fail=1; }
    done
    n=$(wc -l < "$out/m0.f$frame")
    bad=$(diff "$out/m0.f$frame" "$out/m1.f$frame" | grep -c '^<')
    if [ "$n" -gt 0 ] && [ "$bad" = 0 ] && cmp -s "$out/m0.f$frame" "$out/m1.f$frame"; then
        echo "PASS $prot x $dna, $frame frames: $n hits identical"
    else
        echo "FAIL $prot x $dna, $frame frames: $bad of $n hits differ (read protein score i j reverse)"
        diff "$out/m0.f$frame" "$out/m1.f$frame" | head -20
        fail=1
    fi
done
exit "$fail"
