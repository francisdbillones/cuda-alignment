#!/usr/bin/env bash
# Installs what this project needs to build and run on x86_64 Debian or Ubuntu, then test-builds it:
# CUDA 13.4.3, CMake 3.31.12, Ninja 1.13.2 and the distro's g++. Run it as root or with sudo.
# Afterwards, run `source /opt/cuda-alignment/env.sh` in each shell that builds the project.
# It doesn't install the NVIDIA driver, which is needed to run the program.
set -euo pipefail

PREFIX=/opt/cuda-alignment
CUDA_VERSION=13.4.3
CUDA_MD5=d79b49d3ab39f5d45d48d91b31275de4
CMAKE_VERSION=3.31.12
CMAKE_SHA256=0dc2e9a6860f06bf10bd8fadc03e35d9eeb4df46e33763a7e480e987758f385c
NINJA_VERSION=1.13.2
NINJA_SHA256=5749cbc4e668273514150a80e387a957f933c6ed3f5f11e03fb30955e2bbead6

[ "$(uname -sm)" = "Linux x86_64" ] || { echo "This script only supports x86_64 Debian or Ubuntu." >&2; exit 1; }
[ "$(id -u)" -eq 0 ] || { echo "Run this script as root or with sudo." >&2; exit 1; }
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"

apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    g++ curl ca-certificates tar unzip perl libxml2 tmux

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
fetch() { curl -fL --retry 5 --connect-timeout 30 -o "$1" "$2"; }

mkdir -p "$PREFIX"

fetch "$tmp/cuda.run" "https://developer.download.nvidia.com/compute/cuda/$CUDA_VERSION/local_installers/cuda_${CUDA_VERSION}_linux.run"
echo "$CUDA_MD5  $tmp/cuda.run" | md5sum --check -
# Without --silent the installer opens an interactive menu, so follow its log instead.
echo "Installing CUDA $CUDA_VERSION. This can take several minutes."
tail -n 0 -F /var/log/cuda-installer.log 2>/dev/null &
log_tail=$!
installer_status=0
sh "$tmp/cuda.run" --silent --toolkit --toolkitpath="$PREFIX/cuda" --no-man-page --tmpdir="$tmp" || installer_status=$?
kill "$log_tail"
[ "$installer_status" -eq 0 ] || { echo "The CUDA installer failed; see /var/log/cuda-installer.log." >&2; exit "$installer_status"; }
rm "$tmp/cuda.run"

fetch "$tmp/cmake.tgz" "https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-linux-x86_64.tar.gz"
echo "$CMAKE_SHA256  $tmp/cmake.tgz" | sha256sum --check -
mkdir -p "$PREFIX/cmake"
tar -xzf "$tmp/cmake.tgz" --strip-components=1 -C "$PREFIX/cmake"

fetch "$tmp/ninja.zip" "https://github.com/ninja-build/ninja/releases/download/v$NINJA_VERSION/ninja-linux.zip"
echo "$NINJA_SHA256  $tmp/ninja.zip" | sha256sum --check -
unzip -qo "$tmp/ninja.zip" -d "$PREFIX/cmake/bin"

cat >"$PREFIX/env.sh" <<EOF
export PATH="$PREFIX/cuda/bin:$PREFIX/cmake/bin:\$PATH"
export LD_LIBRARY_PATH="$PREFIX/cuda/lib64\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
EOF
. "$PREFIX/env.sh"

cmake -S "$repo" -B "$tmp/build" -G Ninja
cmake --build "$tmp/build"
echo "Test build succeeded."

if nvidia-smi -L 2>/dev/null; then
    echo "NVIDIA driver found. CUDA 13.4.3 programs may need a recent driver; check 'nvidia-smi' if the program won't start."
else
    echo "No working NVIDIA driver found, so the program can build here but not run."
fi
[ -f "$repo/data/BLOSUM62.txt" ] || echo "data/ is missing. It isn't in git, so copy it here before running the program."
