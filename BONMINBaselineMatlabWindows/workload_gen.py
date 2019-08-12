# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 11:03:02 2019

@author: amaity
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def workloadGen(numPhases,numDeadlines,numWorkloads,suffix,deadlineFixed=False):
    """
        Create a workload with "numPhases"
        with different deadlines and workload
        mix and export it in a csv file.
    """
    fl2       = open("workloads-"+suffix+"/wkld_"+str(numPhases)+".csv",'w')
    
    # Write the CSV header
    # fl2.write("Deadline,")
    # for p in range(1,numPhases):
    #     fl2.write("bench-"+str(p)+",")
    # fl2.write("bench-"+str(numPhases)+"\n")
    
    # For Each deadline
    if deadlineFixed :
        deadlines = [14972*numPhases/7]
    else:
        deadlines = np.linspace(600*numPhases,1897*numPhases,numDeadlines)

    for i in range(0,numWorkloads):
        # Create a Random Workloads
        benchid = np.random.randint(5,11,numPhases)

        # Export the workload for each possible deadline
        for d in deadlines:
            # Write the CSV Entry
            fl2.write(str(d)+",")
            for p in range(0,len(benchid)-1):
                fl2.write(str(benchid[p])+",")
            p = p + 1
            fl2.write(str(benchid[p])+"\n")
    
    fl2.close()    

def workloadGenDual(numPhases,numPowers,numWorkloads,suffix,powerFixed=False):
    """
        Create a workload with "numPhases"
        with different (peak) powers and workload
        mix and export it in a csv file.
    """
    fl2       = open("workloads-"+suffix+"/wkld_"+str(numPhases)+".csv",'w')

    
    # For Each deadline
    if powerFixed :
        powers = 12
    else:
        powers = np.linspace(10,30,numPowers)

    for i in range(0,numWorkloads):
        # Create a Random Workloads
        benchid = np.random.randint(5,11,numPhases)

        # Export the workload for each possible deadline
        for d in powers:
            # Write the CSV Entry
            fl2.write(str(d)+",")
            for p in range(0,len(benchid)-1):
                fl2.write(str(benchid[p])+",")
            p = p + 1
            fl2.write(str(benchid[p])+"\n")
    
    fl2.close()   


def exp1():
    """
        Average peak power for applications under PkMin with
        increasing number of phases normalized against their peak
        power under Oracle.
    """
    phases = range(2,51)
    for numPhases in phases:
        workloadGen(numPhases,10,10,"exp1")

def exp3():
    """
        Peak power under PkMin for yyy randomized 20-phase
        real-time applications normalized against their peak power
        under Oracle. Only xxx applications out of hundred have peak
        power higher under PkMin than Oracle, deadline if Fixed
    """
    phases = 20
    workloadGen(phases,-1,120,"exp3",True)

def exp4():
    """
        Peak power under PkMin for a 20-phase real-time
        application with different deadlines normalized against its peak
        power under Oracle. Only six deadlines out of hundred have
        peak power higher under PkMin than Oracle.
    """
    phases = 20
    workloadGen(phases,120,1,"exp4",False)

def exp5():
    """
        Average execution time for applications under EtMin
        normalized against their execution time under Oracle with
        different power budgets.
    """
    phases = range(2,31)
    for numPhases in phases:
        workloadGenDual(numPhases,10,10,"exp5")

def exp6():
    """
        Runtime of PkMin vs BONMIN for applications with
        different number of phases.
    """
    phases = range(31,101)
    for numPhases in phases:
        workloadGen(numPhases,10,10,"exp6",False)

if __name__=="__main__":
    # workloadGen(6,10,10)
    # os.system('mkdir workloads-exp3 workloads-exp4')
    # exp3()
    # exp4()
    os.system('mkdir workloads-exp6')
    # exp5()
    exp6()