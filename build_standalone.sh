#!/bin/sh

HERE=$(pwd)
BUMPS_DIR=$HOME/dev/bumps
REFL1D_DIR=$HOME/dev/refl1d
MOLGROUPS_DIR=$HOME/dev/molgroups
TARGET=$HERE/public/wheels
export BUILD_EXTENSION=True

cd $BUMPS_DIR && pyodide build
#cd $REFL1D_DIR && python setup.py bdist_wheel
cd $REFL1D_DIR && pyodide build
cd $MOLGROUPS_DIR && python setup.py bdist_wheel

cd $HERE
mkdir -p $TARGET
cp $BUMPS_DIR/dist/*.whl $TARGET
cp $REFL1D_DIR/dist/*.whl $TARGET
cp $MOLGROUPS_DIR/dist/*.whl $TARGET

npm run build -- --sourcemap