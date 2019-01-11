#!/bin/bash

exe="${HOME}/Desktop/test2"
reset && g++ -O3 -pg -gdwarf-3 -std=c++11 -Wall ptss_gsltest.cpp -o ${exe} -lm
#reset && g++ -std=c++11 -Wall tmp.cpp -o test2 -lm
#./${exe}
