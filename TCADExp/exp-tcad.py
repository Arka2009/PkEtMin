#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 17:39:42 2019

@author: amaity)
"""
import os

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


if __name__=="__main__":
    write_build(4,False)