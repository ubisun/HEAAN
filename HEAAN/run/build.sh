#!/bin/bash

cd ../lib
make -j4 all
cd ../run
make
