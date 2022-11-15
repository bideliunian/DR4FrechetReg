#!/bin/bash

for args in `seq 1 100`;
do
qsub dist_model12.PBS -v "args=$args"
done
