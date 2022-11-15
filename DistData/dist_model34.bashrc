#!/bin/bash

for args in `seq 1 100`;
do
qsub dist_model34.PBS -v "args=$args"
done
