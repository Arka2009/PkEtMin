# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 16:40:48 2019

@author: amaity
"""
import numpy as np
import pandas as pd

df   = pd.read_csv("profile-lacedfs-sniper_000.csv")
PH1  = df['alloc'].values
PH2  = df['alloc'].values
muet = df['mu-time'].values
ppgp = df['ppg-power'].values
Z    = open("et-power2.csv","w")

for x in PH1:
    for y in PH2:
#        z = (muet[x-1]*ppgp[x-1] + 33*muet[y-1]*ppgp[y-1])/(muet[x-1]+muet[y-1])
        z = np.max([ppgp[x-1],ppgp[y-1]])
        Z.write(f'{z},')
    Z.write("\n")
Z.close()