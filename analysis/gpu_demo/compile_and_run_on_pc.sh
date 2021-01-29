#!/bin/bash
PGFORTRAN=/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/bin/pgfortran
${PGFORTRAN} -acc -mp=nonuma -O3 -ta=tesla:cc60 -Mcuda -Minfo=accel gpu_demo.F90 -o gpu_demo && ./gpu_demo

