#!/bin/bash

make

mkdir -p "${PREFIX}/bin"
cp reads_sorting "${PREFIX}/bin/"
cp sort_the_reads.py "${PREFIX}/bin/"
