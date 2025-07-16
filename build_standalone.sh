#!/bin/sh

HERE=$(pwd)
BUMPS_DIR=$HERE/bumps
REFL1D_DIR=$HERE/refl1d
MOLGROUPS_DIR=$HERE/molgroups
TARGET=$HERE/public/wheels
export BUILD_EXTENSION=True

cd $BUMPS_DIR && rm -rf dist && pyodide build
#cd $REFL1D_DIR && python setup.py bdist_wheel
cd $REFL1D_DIR && rm -rf dist && pyodide build
cd $MOLGROUPS_DIR && rm -rf dist && pyodide build

cd $HERE
mkdir -p $TARGET
rm $TARGET/*.whl
cp $BUMPS_DIR/dist/*.whl $TARGET
cp $REFL1D_DIR/dist/*.whl $TARGET
cp $MOLGROUPS_DIR/dist/*.whl $TARGET

npm run build -- --sourcemap
