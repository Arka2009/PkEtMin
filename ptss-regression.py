#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 22:55:34 2019

@author: amaity
"""

import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt

def main(dir2,bench):
    fl2 = dir2 + "/lace-characteristics-"+bench+"-003.csv"
    df2 = pd.read_csv(fl2)
    a   = df2['alloc'].values
    a2  = [np.log(u) for u in a]
    e   = [1/u for u in df2['latency(M)'].values]
    e2  = df2['latency(M)'].values
    p   = df2['PKG'].values
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(a,e)
    print(f'Exec Time - Bench : {bench}, slope : {slope}, intercept : {intercept}, r : {r_value}')
    slope, intercept, r_value, p_value, std_err = stats.linregress(a,p)
    print(f'Power - Bench : {bench}, slope : {slope}, intercept : {intercept}, r : {r_value}')

    # Plot
    plt.plot(a,e)
    plt.xlabel("Allocation")
    plt.ylabel("Inverse of Execution Time")
    plt.title("Execution Time Characteristics")
    fl3 = dir2+"/lace-lat-viz-"+bench+"-003.pdf"
    plt.savefig(fl3)
    plt.close()

    plt.plot(a,p)
    plt.xlabel("Allocation")
    plt.ylabel("PKP Power in Watts")
    plt.title("Power Characteristics")
    fl3 = dir2+"/lace-power-viz-"+bench+"-003.pdf"
    plt.savefig(fl3)
    plt.close()

    
    
    print("\n")

if __name__=="__main__":
    main("profile-lace-003","cilksort")
    main("profile-lace-003","dfs")
    main("profile-lace-003","fib")
    main("profile-lace-003","pi")
    main("profile-lace-003","queens")