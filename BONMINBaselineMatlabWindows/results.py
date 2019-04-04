# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 07:48:17 2019

@author: amaity
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def exp1():
    """
        Average peak power for applications under PkMin with
        increasing number of phases normalized against their peak
        power under Oracle.
    """
    phases = range(2,31)
    fil3 = open("exp1_bonmin.csv",'w')
    fil3.write("numPhases,ratio,ratio2,etBONMIN,etPkMin\n")
    for numPhases in phases:
        colNames = ["alloc-"+str(u) for u in range(1,numPhases+1)]
        colNames.append("pkp")
        colNames.append("status")
        colNames.append("OH")
        
        fld1 = open("workloads-exp1/wkld_"+str(numPhases)+"_matlab.out.csv","r")
        df1  = pd.read_csv(fld1,header=None,names=colNames)
        df11 = df1[df1.status == "passed"]
        arr1 = df11['pkp'].values

        colNames.append("wrst")
        fld2 = open("workloads-exp1/wkld_"+str(numPhases)+"_cpp.out.csv","r")
        df2  = pd.read_csv(fld2,header=None,names=colNames)
        df21 = df2[df2.status == "passed"]
        arr2 = df21['pkp'].values
        arr3 = df21['wrst'].values
        
        ratio  = np.divide(arr2,arr1)
        ratio2 = np.divide(arr3,arr2)

        etbonmin = df11['OH'].values
        etPkMin  = df21['OH'].values
        fil3.write(f"{numPhases},{np.mean(ratio)},{np.mean(ratio2)},{np.mean(etbonmin)},{np.mean(etPkMin)}\n")
  
    
    fil3.close()

def exp3():
    """
        Peak power under PkMin for yyy randomized 20-phase
        real-time applications normalized against their peak power
        under Oracle. Only xxx applications out of hundred have peak
        power higher under PkMin than Oracle, deadline if Fixed
    """
    
    fil3 = open("exp3_bonmin.csv",'w')
    fil3.write("WorkloadId,ratio\n")
    numPhases = 20
    
    colNames = ["alloc-"+str(u) for u in range(1,numPhases+1)]
    colNames.append("pkp")
    colNames.append("status")
    colNames.append("OH")
    
    fld1 = open("workloads-exp3/wkld_"+str(numPhases)+"_matlab.out.csv","r")
    df1  = pd.read_csv(fld1,header=None,names=colNames)
    df11 = df1[df1.status == "passed"]
    arr1 = df11['pkp'].values
    # print(arr1)

    fld2 = open("workloads-exp3/wkld_"+str(numPhases)+"_cpp.out.csv","r")
    df2  = pd.read_csv(fld2,header=None,names=colNames)
    df21 = df2[df2.status == "passed"]
    arr2 = df21['pkp'].values
    # print(arr2)
    
    ratio = np.divide(arr2,arr1)
    for i in range(1,len(ratio)+1):
        fil3.write(f"{i},{ratio[i-1]}\n")
  
    fil3.close()


def exp4():
    """
        Peak power under PkMin for a 20-phase real-time
        application with different deadlines normalized against its peak
        power under Oracle. Only six deadlines out of hundred have
        peak power higher under PkMin than Oracle.
    """
    
    fil3 = open("exp4_bonmin.csv",'w')
    fil3.write("WorkloadId,ratio\n")
    numPhases = 20
    
    colNames = ["alloc-"+str(u) for u in range(1,numPhases+1)]
    colNames2 = ["deadlines"] + colNames
    colNames.append("pkp")
    colNames.append("status")
    colNames.append("OH")

    fld0 = open("workloads-exp4/wkld_"+str(numPhases)+".csv","r")
    df0  =  pd.read_csv(fld0,header=None,names=colNames2)   
    deadlines = df0['deadlines'].values

    fld1 = open("workloads-exp4/wkld_"+str(numPhases)+"_matlab.out.csv","r")
    df1  = pd.read_csv(fld1,header=None,names=colNames)
    df11 = df1[df1.status == "passed"]
    arr1 = df11['pkp'].values
    # print(arr1)

    fld2 = open("workloads-exp4/wkld_"+str(numPhases)+"_cpp.out.csv","r")
    df2  = pd.read_csv(fld2,header=None,names=colNames)
    df21 = df2[df2.status == "passed"]
    arr2 = df21['pkp'].values
    # print(arr2)
    
    ratio = np.divide(arr2,arr1)
    for i in range(1,len(ratio)+1):
        fil3.write(f"{deadlines[i]},{ratio[i-1]}\n")
  
    fil3.close()

def main():
    exp1()
    # exp3()
    # exp4()

if __name__=="__main__":
    main()