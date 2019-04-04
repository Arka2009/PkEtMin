#!/bin/bash


for i in $(seq 0 4); do 
	grep -w "ufhew4r4" dump3/exp3-deadline${i}.log >> dump3/exp3-pkp.csv; 
done
