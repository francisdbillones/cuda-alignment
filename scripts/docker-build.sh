#!/usr/bin/env bash
# Builds the project inside Docker and writes the output to build/ on the host.
# Extra arguments go to `cmake --build`, e.g. scripts/docker-build.sh --clean-first
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
dockerfile="$repo/docker/Dockerfile"

# The tag follows the Dockerfile's contents, so editing it triggers a rebuild.
image="cuda-alignment-dev:$(shasum -a 256 "$dockerfile" | cut -c1-12)"
if ! docker image inspect "$image" >/dev/null 2>&1; then
    docker build -t "$image" -f "$dockerfile" "$repo/docker"
fi

# The repo is mounted at its macOS path so compiler errors and the CMake cache
# use paths that are also valid on the host.
docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$repo:$repo" \
    --workdir "$repo" \
    "$image" \
    bash -c 'cmake -S . -B build -G Ninja && cmake --build build "$@"' _ "$@"
