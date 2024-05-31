#!/bin/sh

HERE=$(pwd)
BUMPS_DIR=/Users/bbm/dev/bumps
REFL1D_DIR=/Users/bbm/dev/refl1d
TARGET=$HERE/public/wheels
export BUILD_EXTENSION=True

cd $BUMPS_DIR && python setup.py bdist_wheel
export BUMPS_WHEEL_FILE=$(ls dist)
#cd $REFL1D_DIR && python setup.py bdist_wheel
cd $REFL1D_DIR && pyodide build
export REFL1D_WHEEL_FILE=$(ls dist)

cd $HERE
npm run build -- --sourcemap
mkdir -p $TARGET
cp $BUMPS_DIR/dist/*.whl $TARGET
cp $REFL1D_DIR/dist/*.whl $TARGET
