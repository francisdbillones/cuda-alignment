#!/usr/bin/env bash
set -uo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
bin="${1:-$repo/build/cuda_alignment}"
golden="$repo/tests/golden"
out="$(mktemp -d)"
trap 'rm -rf "$out"' EXIT

[ -x "$bin" ] || { echo "No executable at $bin" >&2; exit 2; }

fail=0
for ds in 'testProt testDNA' 'ffsub_prot ffsub_dna'; do
    set -- $ds
    for mode in 0 1; do for frame in 3 6; do for top in 0 1; do
        tag="$1.m$mode.f$frame.t$top"
        printf '%s\n' "$1" "$2" "$mode" "$frame" "$top" | timeout 900 "$bin" > "$out/$tag.raw" 2>&1
        rc=$?
        grep -vE 'Time in ms|Total runtime' "$out/$tag.raw" > "$out/$tag.txt"
        if [ "$rc" = 0 ] && cmp -s "$golden/$tag.txt" "$out/$tag.txt"; then
            echo "PASS $tag"
        else
            echo "FAIL $tag (exit $rc)"
            diff "$golden/$tag.txt" "$out/$tag.txt" | head -20
            fail=1
        fi
    done; done; done
done

echo
if [ "$fail" = 0 ]; then echo 'ALL PASS'; else echo 'SOME FAILED'; fi
exit "$fail"
