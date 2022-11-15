#!/bin/bash

for args in `seq 1 99`;
do
qsub sphere.PBS -v "args=$args"
done
