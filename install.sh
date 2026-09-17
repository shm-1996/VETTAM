#!/usr/bin/env bash
# Install the VETTAM units into a FLASH source tree.
#
#   ./install.sh /path/to/flash/source           copy the units into place
#   ./install.sh --link /path/to/flash/source    symlink them instead
#
# Only the two unit directories are created or replaced. No FLASH core file is
# touched.
set -euo pipefail

MODE=copy
if [ "${1:-}" = "--link" ]; then
    MODE=link
    shift
fi

if [ $# -ne 1 ]; then
    echo "usage: $0 [--link] /path/to/flash/source" >&2
    exit 2
fi

TARGET=$1
REPO=$(cd "$(dirname "$0")" && pwd)

if [ ! -d "$TARGET" ]; then
    echo "error: $TARGET is not a directory" >&2
    exit 1
fi
if [ ! -d "$TARGET/physics/RadTrans/RadTransMain" ]; then
    echo "error: $TARGET does not look like a FLASH source/ directory" >&2
    echo "       (expected physics/RadTrans/RadTransMain beneath it)" >&2
    exit 1
fi

UNITS="physics/RadTrans/RadTransMain/VETTAM physics/materialProperties/Opacity/RadTrans"

for unit in $UNITS; do
    dest="$TARGET/$unit"
    parent=$(dirname "$dest")
    mkdir -p "$parent"
    if [ -e "$dest" ] || [ -L "$dest" ]; then
        backup="$dest.bak.$(date +%Y%m%d%H%M%S)"
        echo "  existing $unit -> $(basename "$backup")"
        mv "$dest" "$backup"
    fi
    if [ "$MODE" = link ]; then
        ln -s "$REPO/$unit" "$dest"
        echo "  linked  $unit"
    else
        cp -R "$REPO/$unit" "$dest"
        echo "  copied  $unit"
    fi
done

echo
echo "Done. Next:"
echo "  1. Point your site Makefile.h at PETSc (PETSC_PATH, FFLAGS_PETSC, LIB_PETSC)."
echo "  2. Add to your simulation Config:"
echo "       REQUIRES physics/RadTrans/RadTransMain/VETTAM"
echo "  3. ./setup YourSimulation -auto ... && cd objdir && make"
