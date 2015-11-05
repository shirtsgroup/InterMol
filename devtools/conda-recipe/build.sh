#!/bin/bash

cp -r $RECIPE_DIR/../.. $SRC_DIR
$PYTHON setup.py clean
$PYTHON setup.py install --single-version-externally-managed --record=record.txt
