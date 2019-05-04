#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import heapq
import queue
import timeit
import pandas as pd


def main():
    fl2  = "build/DistTest.csv"

    df2  = pd.read_csv(fl2)
    
    for i in range(0,8):
        srs  = df2[f'Alloc-{i}'].values
        plt.hist(srs,bins=1000,density=True)
        
    plt.savefig('tmp.pdf')
    plt.xlabel('Execution Time')
    plt.ylabel('Probability Density')
    plt.title('Execution Time Variation of 2D Convolution')
    

if __name__=="__main__":
    main()