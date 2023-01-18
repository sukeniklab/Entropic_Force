#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 16:24:47 2021

This script is an example of how I conduct ensemble analysis for the paper.

I developed a python library to help me analyze multi-trajectory ensembles
which may be shared in my personal Github.

The core function used here is the same as that library.

@author: lemoncatboy
"""
import mdtraj as md
import os, sys
import statistics as st
import numpy as np
import pandas as pd
import re

# Calculate the Hydrogen-bond number using Mdtraj library
# Here the r is the trajectory that only containing the protein atoms
# j is the number of frames in the trajectory.
def calc_HB(r,j,length):
    #Calculate the Hydogen_bond
    hbo=md.wernet_nilsson(r)
    # This is a couter for the loop over frames
    op=0
    # A list to store the final result.
    HB=[]
    # looping over each frame
    while op<j:
        # hbo[op] is a sub numpy array stored donor and acceptor 
        # information of the H-bonds. The length of the array is the number of
        # H-bonds
        if (len(hbo[op])):
            HB.append(len(hbo[op])/length)
        op+=1
    # Average H-bonds per residue per frame
    return(sum(HB)/j)

# Calculate the helical propensity fo the protein.
# Here the t is the full trajectory file(We can try to use r instead)
def calc_Heli(t,length):
    # Calculate the secondary structure using dssp
    dssp = md.compute_dssp(t)
    # Create an empty array with the row number equals to the number of residues
    # We may not need the number 1 here
    dssp_count = np.zeros((1, t.n_residues))
    # Loop over frames
    for i in range(t.n_frames):
        # Loop over residues
        for j in range(t.n_residues):
            # In dssp, letter 'H' means alpha-helix
            if dssp[i,j] in 'H':
                dssp_count[0,j] += 1
    # Get the average of the dssp_count
    dssp_prob = np.divide(dssp_count,t.n_frames)
    summary=dssp_prob.sum(axis=1)
    return(summary[0]/length)

#calc_beta is added 03/28/2022
# Same algorithm as alpha
def calc_Beta(t,length):
    dssp = md.compute_dssp(t)
    dssp_count = np.zeros((1, t.n_residues))
    for i in range(t.n_frames):
        for j in range(t.n_residues):
            if dssp[i,j] in 'E':
                dssp_count[0,j] += 1
    dssp_prob = np.divide(dssp_count,t.n_frames)
    summary = dssp_prob.sum(axis=1)
    return(summary[0]/length)

# Average radius of gyration of the ensemble
def calc_Rg(r):
    d = md.compute_rg(r)
    return(st.mean(d))
# Average end to end distance of the ensemble
# The defination of the end to end distnace is the distance between 
# the first alpha carbon and the last alpha carbon of the protein sequence.

def calc_Re(r,t):
    topology=t.topology
    rpology=topology.select_atom_indices(selection='alpha')
    d = md.compute_distances(r,[[rpology[0],rpology[-1]]])
    listtemp=[]
    for temp in d:
        listtemp.append(float(temp[0]))
    return(st.mean(listtemp))
 
# Get the average of five trajectories. 
# Here, we calculate the weighted average and standard deviation by
# frame numbers     
def average_function(mean_value,stlist,frametraj_temp,meantraj_temp,number_of_repeats):
    framesum=0
    meanadj=0
    for nn in range(0,number_of_repeats,1):
        framesum+=frametraj_temp[nn]
        meanadj+=frametraj_temp[nn]*meantraj_temp[nn]
    meanref=meanadj/framesum
    mean_value.append(meanref)
    stdadj=0
    for nn in range(0,number_of_repeats,1):
        stdadj+=frametraj_temp[nn]*(meantraj_temp[nn]-meanref)*(meantraj_temp[nn]-meanref)
    stlist.append(np.sqrt((stdadj/framesum)*number_of_repeats/(number_of_repeats-1)))

#Read sequence from the fast file
def read_seq():
    seqopen=open('seq.fasta','r')
    seq=seqopen.read()
    while (re.search('\s',seq[-1]) != None):
        seq=seq[:-1]
    x=len(seq)
    return seq,x

# A protocol function for calculating all quantities above.
def easy_standard(pwd,pdb_name,traj_name=False):
    os.chdir(pwd)
    print (pwd)
    seq,length=read_seq()
    if traj_name:
        t = md.load(traj_name,top=pdb_name)
    else:
        t = md.load(pdb_name)
    # Select the protein in the simulation trajectory
    u_slice=t.top.select('protein')
    # Slice the trajectory
    r_sliced=t.atom_slice(u_slice)
    j = r_sliced.n_frames
    return (calc_HB(r_sliced,j,length),
            calc_Heli(t,length),
            calc_Beta(t,length),
            calc_Rg(r_sliced),
            calc_Re(r_sliced,t))
        

if __name__=="__main__":
    pwd=r'F:\DATA_F\Entropy_Force_raw_data\Analysis_script\Example\1_pdb+traj'
    (HB_result,
     Heli_result,
     Beta_result,
     Rg_result,
     Re_result) = easy_standard(pwd,
                                pdb_name='GS16.pdb',
                                traj_name='GS16.xtc')
