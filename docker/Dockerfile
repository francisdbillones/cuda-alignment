# CUDA build environment. CMake and Ninja match the versions on the macOS host.
FROM nvidia/cuda:13.4.3-devel-ubuntu24.04

ARG CMAKE_VERSION=3.31.12
ARG NINJA_VERSION=1.13.2

RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends ca-certificates curl unzip; \
    rm -rf /var/lib/apt/lists/*; \
    arch="$(uname -m)"; \
    curl -fsSL "https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-${arch}.tar.gz" \
        | tar -xz --strip-components=1 -C /usr/local; \
    case "$arch" in \
        aarch64) ninja_zip=ninja-linux-aarch64.zip ;; \
        *) ninja_zip=ninja-linux.zip ;; \
    esac; \
    curl -fsSL -o /tmp/ninja.zip "https://github.com/ninja-build/ninja/releases/download/v${NINJA_VERSION}/${ninja_zip}"; \
    unzip -q /tmp/ninja.zip -d /usr/local/bin; \
    rm /tmp/ninja.zip

# Skip the image's licence banner and driver warning on every build.
ENTRYPOINT []
