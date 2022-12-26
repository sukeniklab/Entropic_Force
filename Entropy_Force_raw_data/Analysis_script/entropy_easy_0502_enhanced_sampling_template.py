#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 16:24:47 2021

@author: lemoncatboy
"""

import mdtraj as md
import os,sys
import numpy as np
import pandas as pd
import re
import entropy_library as el
import shape_protocol_library as spl

#Read sequence from the fast file
def read_seq():
    seqopen=open('seq.fasta','r')
    seq=seqopen.read()
    while (re.search('\s',seq[-1]) != None):
        seq=seq[:-1]
    x=len(seq)
    return seq,x

# A protocol function for calculating entropy with a surface constratint 
# and enhanced sampling algorithm.
# Compared to other scripts, here the angle_theta = (np.pi-angle_theta)
def enhanced_sampling_surface(pwd,distance_d,angle_theta,pdb_name,traj_name=False):
    os.chdir(pwd)
    print (pwd)
    seq,length=read_seq()
    # Load the trajectory
    if traj_name:
        t = md.load(traj_name,top=pdb_name)
    else:
        t = md.load(pdb_name)
    u_slice=t.top.select('protein')
    r_sliced=t.atom_slice(u_slice)
    topology=t.topology
    r_alpha=topology.select_atom_indices(selection='alpha')
    # The angle defination is different so we need to use np.pi/2-angle_theta
    result=spl.compute_forbiden_rotation_enhanced(distance_d,np.pi/2-angle_theta,t,r_alpha,t.n_frames)
    np.savetxt('Entropy'+str(round(np.pi/2-angle_theta,2))+'.csv',result, delimiter=',')
    ratio=np.average(result)   
    return ratio, t.n_frames

if __name__=="__main__":
    #pwd=r'F:\DATA_F\Entropy_Force_raw_data\Analysis_script\Example\1_pdb+traj'
    pwd=r'F:\DATA_F\Entropy_Force_raw_data\Analysis_script\Example\2_pdb_only'
    distance_d=0
    angle_theta=np.pi / 2 - np.pi / 180
    ratio,frame=enhanced_sampling_surface(pwd,distance_d,angle_theta,pdb_name='GS16.pdb')
    print("$\Delta S/k$ = {:.2f}, {:.2f} frames".format(ratio,frame))