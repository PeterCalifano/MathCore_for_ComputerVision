FROM mcr.microsoft.com/devcontainers/cpp:1-ubuntu-24.04

ARG REINSTALL_CMAKE_VERSION_FROM_SOURCE="none"

# Optionally replace the distribution CMake with a pinned source build.
COPY ./reinstall-cmake.sh /tmp/reinstall-cmake.sh
RUN if [ "${REINSTALL_CMAKE_VERSION_FROM_SOURCE}" != "none" ]; then \
        chmod +x /tmp/reinstall-cmake.sh \
        && /tmp/reinstall-cmake.sh "${REINSTALL_CMAKE_VERSION_FROM_SOURCE}"; \
    fi \
    && rm -f /tmp/reinstall-cmake.sh

# Install MathCore development dependencies. Valgrind is optional because it
# is a profiling aid, not a build requirement.
ARG INSTALL_VALGRIND="on"
ARG VALGRIND_VERSION="3.27.1"
COPY ./custom-setup.sh /tmp/custom-setup.sh
RUN chmod +x /tmp/custom-setup.sh \
    && INSTALL_VALGRIND="${INSTALL_VALGRIND}" VALGRIND_VERSION="${VALGRIND_VERSION}" /tmp/custom-setup.sh \
    && rm -f /tmp/custom-setup.sh

# Standalone image builds can request CUDA from NVIDIA's apt repository. The
# devcontainer flow normally supplies CUDA through the nvidia-cuda feature.
ARG INSTALL_CUDA="off"
ARG CUDA_VERSION="12.9"
COPY ./cuda-setup.sh /tmp/cuda-setup.sh
RUN chmod +x /tmp/cuda-setup.sh \
    && INSTALL_CUDA="${INSTALL_CUDA}" CUDA_VERSION="${CUDA_VERSION}" /tmp/cuda-setup.sh \
    && rm -f /tmp/cuda-setup.sh

ENV CUDA_HOME=/usr/local/cuda
ENV PATH=/usr/local/cuda/bin:${PATH}
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

# Make GPG signing work in interactive shells without evaluating tty during
# the non-interactive image build.
RUN echo "export GPG_TTY=\$(tty)" >> /etc/bash.bashrc
