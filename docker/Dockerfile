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

# clang-tidy 23.1.4 (same LLVM as the macOS host), from the statically linked PyPI wheel.
# It bundles clang's resource headers, including the CUDA runtime wrapper, and avoids
# installing LLVM from the slow ports.ubuntu.com mirror on arm64.
RUN set -eux; \
    case "$(uname -m)" in \
        aarch64) \
            url=https://files.pythonhosted.org/packages/7d/95/ba9463f94319dddd0c3a9d3f65d432ddf71fe9b203dcd7ea63bcff8999e2/clang_tidy-23.1.4-py2.py3-none-manylinux_2_26_aarch64.manylinux_2_28_aarch64.whl; \
            sha=87596cd77dedfdea23ced1d13e57e869dcaee99745bfe942796196da3c197671 ;; \
        x86_64) \
            url=https://files.pythonhosted.org/packages/80/85/9b483f300124a9245ec6da8975612f1841e26da13b27ad8d17f2006e13bd/clang_tidy-23.1.4-py2.py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl; \
            sha=6bb31e4f319ecbe14acb75f23c3050d0f1d4b1d92ea1c66091bc76e4eac2c1ac ;; \
        *) echo "no clang-tidy wheel for $(uname -m)" >&2; exit 1 ;; \
    esac; \
    curl -fsSL -o /tmp/clang-tidy.whl "$url"; \
    echo "$sha  /tmp/clang-tidy.whl" | sha256sum -c -; \
    unzip -q /tmp/clang-tidy.whl 'clang_tidy/data/*' -d /tmp/clang-tidy; \
    mkdir -p /opt/clang-tidy; \
    cp -r /tmp/clang-tidy/clang_tidy/data/. /opt/clang-tidy/; \
    ln -s /opt/clang-tidy/bin/clang-tidy /usr/local/bin/clang-tidy; \
    rm -rf /tmp/clang-tidy /tmp/clang-tidy.whl

RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends git; \
    rm -rf /var/lib/apt/lists/*

# Skip the image's licence banner and driver warning on every build.
ENTRYPOINT []
