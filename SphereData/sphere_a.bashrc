#!/bin/bash

for args in `seq 1 100`;
do
qsub sphere_a.PBS -v "args=$args"
done
