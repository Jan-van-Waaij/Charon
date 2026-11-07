#!/usr/bin/env bash
# Build, test and package a release of the static charon executable.
#
#   cpp/build-release.sh            # uses the system Eigen and zlib
#   cpp/build-release.sh --fetch    # downloads Eigen and zlib (no system packages needed)
#
# The package (charon-<version>-<os>-<arch>.tar.gz) and its SHA-256 checksum are
# written to dist/ in the repository root. On macOS the executable cannot be fully
# static; there zlib is downloaded and linked statically.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
build="$here/build-release"
dist="$root/dist"

fetch=OFF
static=ON
if [[ "${1:-}" == "--fetch" ]]; then
    fetch=ON
fi
if [[ "$(uname -s)" == "Darwin" ]]; then
    # macOS has no fully static executables; link zlib statically instead.
    static=OFF
    fetch=ON
fi

jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"

cmake -S "$here" -B "$build" -DCMAKE_BUILD_TYPE=Release -DCHARON_STATIC="$static" -DCHARON_FETCH_DEPS="$fetch"
cmake --build "$build" --parallel "$jobs"
(cd "$build" && ctest --output-on-failure)

if command -v ldd >/dev/null 2>&1 && ldd "$build/charon" >/dev/null 2>&1; then
    echo "warning: $build/charon is dynamically linked:" >&2
    ldd "$build/charon" >&2
fi

mkdir -p "$dist"
(cd "$build" && cpack -B "$dist" >/dev/null)
rm -rf "$dist/_CPack_Packages"
cd "$dist"
for package in charon-*.tar.gz charon-*.zip; do
    [[ -e "$package" ]] || continue
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$package" > "$package.sha256"
    else
        shasum -a 256 "$package" > "$package.sha256"
    fi
    echo "Created $dist/$package"
done
"$build/charon" --version
