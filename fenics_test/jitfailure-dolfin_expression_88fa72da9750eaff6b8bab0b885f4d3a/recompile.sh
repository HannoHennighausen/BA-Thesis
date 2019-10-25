#!/bin/bash
# Execute this file to recompile locally
/home/bq_hhennighausen/anaconda3/envs/fenicsproject/bin/x86_64-conda_cos6-linux-gnu-c++ -Wall -shared -fPIC -std=c++11 -O3 -fno-math-errno -fno-trapping-math -ffinite-math-only -I/home/bq_hhennighausen/anaconda3/envs/fenicsproject/include -I/home/bq_hhennighausen/anaconda3/envs/fenicsproject/include/eigen3 -I/home/bq_hhennighausen/anaconda3/envs/fenicsproject/.cache/dijitso/include dolfin_expression_88fa72da9750eaff6b8bab0b885f4d3a.cpp -L/home/bq_hhennighausen/anaconda3/envs/fenicsproject/lib -L/home/bq_hhennighausen/anaconda3/envs/fenicsproject/home/bq_hhennighausen/anaconda3/envs/fenicsproject/lib -L/home/bq_hhennighausen/anaconda3/envs/fenicsproject/.cache/dijitso/lib -Wl,-rpath,/home/bq_hhennighausen/anaconda3/envs/fenicsproject/.cache/dijitso/lib -lmpi -lmpicxx -lpetsc -lslepc -lz -lhdf5 -lboost_timer -ldolfin -olibdijitso-dolfin_expression_88fa72da9750eaff6b8bab0b885f4d3a.so