# Built and pushed by the build_image_sycl CI job (Kaniko) to $SYCL_IMAGE in the
# GitLab Container Registry; that job fires on any change to this file. Bump the
# SYCL_IMAGE tag counter in .gitlab-ci.yml when editing here.
# Pinned hpckit tag: https://hub.docker.com/r/intel/oneapi-hpckit/tags
FROM intel/oneapi-hpckit:2025.3.1-0-devel-ubuntu24.04

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      openmpi-bin libopenmpi-dev \
      libadios2-mpi-c++-dev libadios2-mpi-c++11-dev libadios2-mpi-core-dev \
      libadios2-common-c++11-dev libadios2-common-core-dev \
 && rm -rf /var/lib/apt/lists/*
