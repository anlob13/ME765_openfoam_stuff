#!/bin/bash
set -e
cd "$(dirname "$0")"

export WM_PROJECT_DIR=/usr/lib/openfoam/openfoam2412
export WM_PROJECT=OpenFOAM
export WM_PROJECT_VERSION=2412
export FOAM_ETC=$WM_PROJECT_DIR/etc
export WM_OPTIONS=linux64GccDPInt32Opt
FOAM_LIBBIN=$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib
export LD_LIBRARY_PATH=$FOAM_LIBBIN:$FOAM_LIBBIN/dummy:$LD_LIBRARY_PATH

INCLUDE_DIRS="-I$WM_PROJECT_DIR/src/OpenFOAM/lnInclude -I$WM_PROJECT_DIR/src/OSspecific/POSIX/lnInclude"

g++ -std=c++14 -DWM_DP -DWM_LABEL_SIZE=32 -DNoRepository $INCLUDE_DIRS \
    mlForwardPass_full.C regime_wall_model.C classifier_ensemble.C unified_wall_model.C test_main.C \
    -L"$FOAM_LIBBIN" -L"$FOAM_LIBBIN/dummy" -lOpenFOAM -lPstream \
    -Wl,-rpath,"$FOAM_LIBBIN" -Wl,-rpath,"$FOAM_LIBBIN/dummy" \
    -o testUnifiedWallModel

echo "Build OK - running verification..."
./testUnifiedWallModel
