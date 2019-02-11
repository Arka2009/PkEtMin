#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 22:55:34 2019

@author: amaity
"""

import numpy as np
from scipy import stats
import pandas as pd

def main(fl2):
    df2 = pd.read_csv(fl2)
    a   = df2['alloc'].values
    e   = [1/u for u in df2['latency(M)'].values]
    p   = df2['PKG'].values
    
#    slope, intercept, r_value, p_value, std_err = stats.linregress(a,e)
#    print(f'slope : {slope}, intercept : {intercept}, r : {r_value}')
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(a,p)
    print(f'slope : {slope}, intercept : {intercept}, r : {r_value}')
    
    print("\n")

if __name__=="__main__":
    main("profile-lace-003/profile-lace-dfs-003.csv")
    main("profile-lace-003/profile-lace-fib-003.csv")
    main("profile-lace-003/profile-lace-queens-003.csv")
    main("profile-lace-003/profile-lace-pi-003.csv")
    main("profile-lace-003/profile-lace-cilksort-003.csv")