#!/usr/bin/env bash
# Install the VETTAM units into a FLASH source tree.
#
#   ./install.sh /path/to/flash/source              copy the units into place
#   ./install.sh --link /path/to/flash/source       symlink them instead
#   ./install.sh --no-tests /path/to/flash/source   omit the unit test suite
#
# Only the directories listed in UNITS below are created or replaced, each at
# the same path it occupies in this repository. No FLASH core file is touched.
set -euo pipefail

MODE=copy
WITH_TESTS=yes
while [ $# -gt 0 ]; do
    case "$1" in
        --link)     MODE=link;     shift ;;
        --no-tests) WITH_TESTS=no; shift ;;
        -*)         echo "unknown option: $1" >&2; exit 2 ;;
        *)          break ;;
    esac
done

if [ $# -ne 1 ]; then
    echo "usage: $0 [--link] [--no-tests] /path/to/flash/source" >&2
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

# The unit test suite is a simulation directory rather than a physics unit, but
# it sits at its FLASH path here too, so it installs the same way.
if [ "$WITH_TESTS" = yes ]; then
    UNITS="$UNITS Simulation/SimulationMain/unitTest/VETTAM"
fi

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
if [ "$WITH_TESTS" = yes ]; then
    echo
    echo "Verify the build against the analytic solutions (needs numpy and h5py):"
    echo "  $TARGET/Simulation/SimulationMain/unitTest/VETTAM/run_tests.py"
fi
