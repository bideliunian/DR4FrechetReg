#!/bin/bash

for args in `seq 1 100`;
do
qsub spd.PBS -v "args=$args"
done
