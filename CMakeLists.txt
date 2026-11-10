cmake_minimum_required(VERSION 3.18)

project(cuda_alignment
    DESCRIPTION "GPU-accelerated frameshift-aware DNA-to-protein local alignment"
    LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)

# Timings are meaningless without optimization, so default to Release.
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT is_multi_config AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CUDA_ALIGNMENT_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data" CACHE PATH
    "Directory containing BLOSUM62.txt and the input .fasta files")

add_executable(cuda_alignment src/main.cu src/blosum62.cu src/fasta.cpp)

set_target_properties(cuda_alignment PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(cuda_alignment PRIVATE include)

target_compile_definitions(cuda_alignment PRIVATE
    CUDA_ALIGNMENT_DATA_DIR="${CUDA_ALIGNMENT_DATA_DIR}")
