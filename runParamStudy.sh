#!/bin/bash

. .venv/bin/activate
cd build

./runRe*.sh

python3 parameterStudy.py