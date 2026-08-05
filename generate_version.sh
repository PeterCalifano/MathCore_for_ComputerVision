#!/usr/bin/env bash
# Generate VERSION file without building.
# Fallback chain: git tags --> existing VERSION file --> hardcoded default.

set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION_FILE="${SCRIPT_DIR}/VERSION"

# Default version (last resort)
DEFAULT_MAJOR=0
DEFAULT_MINOR=0
DEFAULT_PATCH=0

# Helpers (matching build_lib.sh style)
info() { echo -e "\e[34m[INFO]\e[0m $*"; }
warn() { echo -e "\e[33m[WARN]\e[0m $*" >&2; }
die() { echo -e "\e[31m[ERROR]\e[0m $*" >&2; exit 1; }

usage() {
    cat <<'EOF'
Usage:
  ./generate_version.sh [options]

Options:
  --sync-ros2     Explicitly synchronize ros2/*/package.xml project metadata.
  --no-sync-ros2  Write VERSION without synchronizing ROS 2 package metadata.
  -h, --help      Show this help.

By default, ROS 2 package metadata is synchronized when the supported overlay
helper is present.
EOF
}

# Automatic mode keeps the standalone generator safe in non-ROS repositories,
# while enabling metadata synchronization when the complete helper is present.
sync_ros2=auto

parse_args() {
    while (($# > 0)); do
        case "$1" in
            --sync-ros2)
                sync_ros2=true
                shift
                ;;
            --no-sync-ros2)
                sync_ros2=false
                shift
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
}

parse_args "$@"

version_major=""
version_minor=""
version_patch=""
full_version=""
source=""

# 1. Try git
if command -v git >/dev/null 2>&1 && git -C "$SCRIPT_DIR" rev-parse --git-dir >/dev/null 2>&1; then
    git_tag=$(git -C "$SCRIPT_DIR" describe --tags --always 2>/dev/null || true)
    git_hash=$(git -C "$SCRIPT_DIR" rev-parse --short=7 HEAD 2>/dev/null || true)

    # Strip leading 'v' and match semver
    clean_tag="${git_tag#v}"
    if [[ "$clean_tag" =~ ^([0-9]+)\.([0-9]+)\.([0-9]+)(-.*)? ]]; then
        version_major="${BASH_REMATCH[1]}"
        version_minor="${BASH_REMATCH[2]}"
        version_patch="${BASH_REMATCH[3]}"
        full_version="${version_major}.${version_minor}.${version_patch}+${git_hash}"
        source="git"
    fi
fi

# 2. Try existing VERSION file
if [[ -z "$source" && -f "$VERSION_FILE" ]]; then
    if grep -qP 'Project version: \d+\.\d+\.\d+' "$VERSION_FILE"; then
        line=$(grep -oP 'Project version: \K\d+\.\d+\.\d+' "$VERSION_FILE")
        IFS='.' read -r version_major version_minor version_patch <<< "$line"

        # Try to read full version line
        full_line=$(grep -oP 'Full version: \K.+' "$VERSION_FILE" 2>/dev/null || true)
        if [[ -n "$full_line" ]]; then
            full_version="$full_line"
        else
            full_version="${version_major}.${version_minor}.${version_patch}"
        fi
        source="VERSION file"
    else
        warn "VERSION file exists but could not be parsed"
    fi
fi

# 3. Fallback to hardcoded defaults
if [[ -z "$source" ]]; then
    version_major=$DEFAULT_MAJOR
    version_minor=$DEFAULT_MINOR
    version_patch=$DEFAULT_PATCH
    full_version="${version_major}.${version_minor}.${version_patch}"
    source="hardcoded default"
    warn "No git tags or VERSION file found. Using default version: ${full_version}"
fi

# Write VERSION file
version_string="${version_major}.${version_minor}.${version_patch}"
{
    echo "Project version: ${version_string}"
    echo "Full version: ${full_version}"
} > "$VERSION_FILE"

info "Version ${full_version} (from ${source}) written to ${VERSION_FILE}"

sync_ros2_package_metadata() {
    if [[ "${sync_ros2}" == false ]]; then
        return
    fi

    local ros2_dir_="${SCRIPT_DIR}/ros2"
    if [[ ! -d "${ros2_dir_}" ]]; then
        if [[ "${sync_ros2}" == true ]]; then
            info "ROS 2 overlay not present; skipping package metadata sync"
        fi
        return
    fi

    if [[ "${source}" == "hardcoded default" ]]; then
        warn "Skipping ROS 2 package metadata sync because the version came from the hardcoded default"
        return
    fi

    if [[ ! "${version_string}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
        warn "Skipping ROS 2 package metadata sync because '${version_string}' is not strict X.Y.Z"
        return
    fi

    local metadata_helper_="${ros2_dir_}/tools/sync_package_metadata.py"
    if [[ ! -f "${metadata_helper_}" ]]; then
        if [[ "${sync_ros2}" == true ]]; then
            warn "ROS 2 package metadata helper is missing; skipping metadata sync"
        fi
        return
    fi
    if ! command -v python3 >/dev/null 2>&1; then
        warn "python3 is unavailable; skipping ROS 2 package metadata sync"
        return
    fi

    python3 "${metadata_helper_}" \
        --project-root "${SCRIPT_DIR}" \
        --ros2-dir "${ros2_dir_}" \
        --version "${version_string}"
}

sync_ros2_package_metadata
