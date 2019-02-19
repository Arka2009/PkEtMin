#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 17:39:42 2019

@author: amaity)
"""
import os
import numpy as np

def write_build(nph,shutdown_oracle):
    # Create Defines File
    fd=open("ptss_config_nph.hpp","w")
    fd.write("#ifndef __PTSS_NPH_CONFIG\n")
    fd.write("#define __PTSS_NPH_CONFIG\n")
    fd.write("#define NPH \t"+str(nph)+"\n")
    if shutdown_oracle :
        fd.write("#define SHUTDOWN_ORACLE\n")
    fd.write("#endif")
    fd.close()

    # Build the executable
    os.system("cd /home/amaity/Dropbox/NUS-Research/ptss_risk_model/ptss-dse/build && make && cd -")

def run_exp1():
    """
        Normalized (pkp) along with increasing 
        number of phases (d)
    """
    for d in range(1,7):
        deadline2 = np.linspace(600*d,1897*d,1000)#np.random.uniform(800*d,5000*d,1000) #[u in 3178*d
        write_build(d,False)
        iter = 0
        for deadline in deadline2:
            print("Exp1:Running with #(phase)="+str(d)+", with deadline="+str(deadline))
            cmd = "/home/amaity/Dropbox/NUS-Research/ptss_risk_model/ptss-dse/build/ptssdse "+str(deadline)+" > dump1/exp1-ph"+str(d)+"-deadline"+str(iter)+".log"
            os.system(cmd)
            iter = iter+1

def run_exp2():
    d = 5
    deadline = 9972
    write_build(d,False)
    print("Exp2:Running")
    cmd = "/home/amaity/Dropbox/NUS-Research/ptss_risk_model/ptss-dse/build/ptssdse2 "+str(deadline)+" > dump2/exp2.log"
    os.system(cmd)

def run_exp3():
    d = 5
    deadline2 = np.linspace(600*d,1897*d,100)#np.random.uniform(800*d,5000*d,1000) #[u in 3178*d
    write_build(d,False)
    iter = 0
    for deadline in deadline2:
        print("Exp3:Running with deadline="+str(deadline))
        cmd = "/home/amaity/Dropbox/NUS-Research/ptss_risk_model/ptss-dse/build/ptssdse "+str(deadline)+" > dump3/exp3-deadline"+str(iter)+".log"
        os.system(cmd)
        iter = iter + 1

def run_exp4():
    for d in range(1,64):
        print("Exp4:Running with #(phase)="+str(d)+", Oracle switched Off")
        write_build(d,True)
        deadline = 1088*d
        cmd = "/home/amaity/Dropbox/NUS-Research/ptss_risk_model/ptss-dse/build/ptssdse "+str(deadline)+" > dump4/exp4-ph"+str(d)+".log"
        os.system(cmd)


def main():
    d = 5
    deadline2 = np.linspace(600*d,1897*d,100)
    fl2 = open("deadline2.csv","w")
    fl2.write("deadline\n")
    for d2 in deadline2:
        fl2.write(str(d2)+"\n")
    fl2.close()

if __name__=="__main__":
    # run_exp1()
#    run_exp2()
#    run_exp3()
    # run_exp4()
     main()