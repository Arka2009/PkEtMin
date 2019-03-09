#!/bin/bash

srcs="ptss_bonmin.cpp ptss_test_bonmin.cpp"

g++ -std=c++11 -I../include -I/home/amaity/Desktop/Bonmin-1.8.7/include/coin -L/home/amaity/Desktop/Bonmin-1.8.7/lib ${srcs} -o test2 -lbonmin -lm -lCoinUtils
