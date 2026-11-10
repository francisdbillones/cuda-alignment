#!/usr/bin/env bash
set -uo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
mode="${1:-1}"
bin="${2:-$repo/build/cuda_alignment}"
golden="$repo/tests/golden/hits"
out="$(mktemp -d)"
trap 'rm -rf "$out"' EXIT

[ -x "$bin" ] || { echo "No executable at $bin" >&2; exit 2; }
mkdir -p "$golden"

fail=0
for ds in 'testProt testDNA' 'ffsub_prot ffsub_dna' 'Fruit_Fly_Protein_Final finalInput' \
    'Mouse_Crest_Protein_Final finalInput'; do
    set -- $ds
    for frame in 3 6; do
        tag="$1.$2.f$frame"
        printf '%s\n' "$1" "$2" "$mode" "$frame" 0 |
            HITS_DUMP="$out/$tag.txt" timeout 1800 "$bin" > "$out/$tag.log" 2>&1
        rc=$?
        if [ "$rc" != 0 ] || [ ! -s "$out/$tag.txt" ]; then
            echo "FAIL $tag (exit $rc)"
            tail -5 "$out/$tag.log"
            fail=1
        elif [ "${UPDATE_GOLDEN:-0}" = 1 ]; then
            cmp -s "$golden/$tag.txt" "$out/$tag.txt" && echo "SAME $tag" || echo "UPDATED $tag"
            cp "$out/$tag.txt" "$golden/$tag.txt"
        elif cmp -s "$golden/$tag.txt" "$out/$tag.txt"; then
            echo "PASS $tag ($(wc -l < "$out/$tag.txt") hits)"
        else
            bad=$(diff "$golden/$tag.txt" "$out/$tag.txt" | grep -c '^<')
            echo "FAIL $tag: $bad of $(wc -l < "$golden/$tag.txt") hits differ (read protein score i j reverse)"
            diff "$golden/$tag.txt" "$out/$tag.txt" | head -20
            fail=1
        fi
    done
done

echo
if [ "$fail" = 0 ]; then echo 'ALL PASS'; else echo 'SOME FAILED'; fi
exit "$fail"
