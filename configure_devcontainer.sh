#!/usr/bin/env bash
# Configure MathCore's base image and optional CUDA devcontainer feature.
set -Eeuo pipefail

readonly ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly DEVCONTAINER_DIR="${ROOT_DIR}/.devcontainer"
readonly DEVCONTAINER_JSON="${DEVCONTAINER_DIR}/devcontainer.json"
readonly DOCKERFILE="${DEVCONTAINER_DIR}/Dockerfile"
readonly JSON_WRITER="${DEVCONTAINER_DIR}/update_devcontainer_json.py"
readonly DEFAULT_CUDA_VERSION="12.9"

cuda_mode="preserve"
cuda_version="${DEFAULT_CUDA_VERSION}"
gpu_runtime="auto"
base=""
base_image=""

usage() {
  cat <<'EOF'
Configure the MathCore development container without building it.

Usage: ./configure_devcontainer.sh [options]

Options:
  --cuda                 Enable the NVIDIA CUDA devcontainer feature.
  --no-cuda              Disable the NVIDIA CUDA devcontainer feature.
  --cuda-version VERSION CUDA feature version (default: 12.9).
  --gpu-runtime RUNTIME  GPU passthrough: auto, docker, or podman.
  --base NAME            ubuntu-24.04, ubuntu-22.04, ubuntu-20.04,
                         debian-12, or debian-11.
  --base-image IMAGE     Set an exact devcontainers/cpp-compatible image.
  -h, --help             Show this help.

Examples:
  ./configure_devcontainer.sh --cuda --gpu-runtime docker
  ./configure_devcontainer.sh --no-cuda --base ubuntu-24.04
  ./configure_devcontainer.sh --base-image ghcr.io/example/cpp:latest
EOF
}

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 1
}

resolve_gpu_runtime() {
  local requested="$1"
  if [[ "${requested}" != "auto" ]]; then
    printf '%s\n' "${requested}"
  elif command -v docker >/dev/null 2>&1; then
    printf 'docker\n'
  elif command -v podman >/dev/null 2>&1; then
    printf 'podman\n'
  else
    printf '[WARN] Neither docker nor podman was detected; writing Docker GPU arguments.\n' >&2
    printf 'docker\n'
  fi
}

while (($# > 0)); do
  case "$1" in
    --cuda)
      cuda_mode="on"
      shift
      ;;
    --no-cuda)
      cuda_mode="off"
      shift
      ;;
    --cuda-version)
      (($# >= 2)) || die '--cuda-version requires a value.'
      cuda_version="$2"
      shift 2
      ;;
    --gpu-runtime)
      (($# >= 2)) || die '--gpu-runtime requires auto, docker, or podman.'
      gpu_runtime="$2"
      shift 2
      ;;
    --base)
      (($# >= 2)) || die '--base requires a value.'
      base="$2"
      shift 2
      ;;
    --base-image)
      (($# >= 2)) || die '--base-image requires a value.'
      base_image="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

[[ -f "${DOCKERFILE}" ]] || die "Dockerfile not found: ${DOCKERFILE}"
[[ -f "${JSON_WRITER}" ]] || die "JSON updater not found: ${JSON_WRITER}"
[[ -z "${base}" || -z "${base_image}" ]] || die '--base and --base-image are mutually exclusive.'

case "${gpu_runtime}" in
  auto|docker|podman) ;;
  *) die "Unsupported GPU runtime: ${gpu_runtime}" ;;
esac

if [[ -n "${base}" ]]; then
  case "${base}" in
    ubuntu-24.04|ubuntu-22.04|ubuntu-20.04|debian-12|debian-11)
      base_image="mcr.microsoft.com/devcontainers/cpp:1-${base}"
      ;;
    *) die "Unsupported base: ${base}" ;;
  esac
fi

if [[ -n "${base_image}" ]]; then
  sed -i "1s|^FROM .*|FROM ${base_image}|" "${DOCKERFILE}"
fi

if [[ "${cuda_mode}" == "preserve" ]]; then
  if { command -v rg >/dev/null 2>&1 &&
       rg -q 'devcontainers/features/nvidia-cuda' "${DEVCONTAINER_JSON}"; } ||
     { ! command -v rg >/dev/null 2>&1 &&
       grep -q 'devcontainers/features/nvidia-cuda' "${DEVCONTAINER_JSON}"; }; then
    cuda_mode="on"
  else
    cuda_mode="off"
  fi
fi

resolved_gpu_runtime="$(resolve_gpu_runtime "${gpu_runtime}")"
temporary_json="$(mktemp "${DEVCONTAINER_DIR}/devcontainer.json.XXXXXX")"
trap 'rm -f -- "${temporary_json}"' EXIT

CUDA="${cuda_mode}" \
CUDA_VERSION="${cuda_version}" \
DEVCONTAINER_GPU_RUNTIME="${resolved_gpu_runtime}" \
DEVCONTAINER_JSON_PATH="${DEVCONTAINER_JSON}" \
python3 "${JSON_WRITER}" >"${temporary_json}"
mv -- "${temporary_json}" "${DEVCONTAINER_JSON}"
trap - EXIT

printf '[INFO] Updated %s\n' "${DEVCONTAINER_JSON}"
printf '[INFO] CUDA: %s (%s), GPU runtime: %s\n' \
  "${cuda_mode}" "${cuda_version}" "${resolved_gpu_runtime}"
printf '[INFO] Container image was not built or run.\n'
