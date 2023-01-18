# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 17:43:42 2022

@author: ShaharGroup-fyu
"""
import mdtraj as md
import numpy as np
import sys
import entropy_library as el

# Enhanced sampling
    # Compared to other scripts, here the angle_theta = (np.pi-angle_theta
def compute_forbiden_rotation_enhanced(distance,angle,trajectory,location_of_alpha,number_of_frames):
    current_frame=0
    # Rotation Interval 
    # This number determined
    rot_angle=60
    # How many rotations should be made
    rot_repeat=int(360/rot_angle)-1
    # Convert degree to rad
    rot_angle_pi=rot_angle/360*2*3.141593
    # Create a np array to save memory processing time
    #This line is for the debugging which will give out the residue that is forbiddened
    #result=np.zeros((number_of_frames,rot_repeat+1,2))
    #This line will only give out the True or False value
    result=np.zeros((number_of_frames,rot_repeat+1))
    while current_frame<number_of_frames:
        objectframe=trajectory.xyz[current_frame,:,:]
        # Shift the object to make a1 as the origin
        object_shifted=el.change_origin(objectframe,objectframe[location_of_alpha[0],:])
        a1a2=objectframe[location_of_alpha[1],:]-objectframe[location_of_alpha[0],:]
        #print(a1a2)
        a1a3=objectframe[location_of_alpha[2],:]-objectframe[location_of_alpha[0],:]
        #print(a1a3)
        a1a2a3_norm_vec=el.unit_normal_vector(a1a2,a1a3)
        #print(a1a2a3_norm_vec)
        # Rotate a1a2 by theta along the norm_vec we can get a1a2p
        # Which is parallel to the constraint surface
        a1a2p=el.vector_rotation(a1a2a3_norm_vec,a1a2,angle)
        #print(a1a2p)
        # Here we created a plane that is parallel to the constrain plane and go through the a1
        a1_r_plane_norm_vec=el.unit_normal_vector(a1a2p,a1a2a3_norm_vec)
        #print(a1_r_plane_norm_vec)
        # The distance between a2 and the constraint plane should be larger than given distance
        sign_vector=np.dot(a1a2,a1_r_plane_norm_vec)
        sign=np.sign(sign_vector)
        #now = datetime.now().time()
        #print("now =", now)
        current_rotation=0
        rot_mat=el.rotation_matrix(rot_angle_pi,a1a2)
        conformation_rotated=object_shifted
        while current_rotation <= rot_repeat:
            limit_ak=[True,0]
            #print(current_rotation)
            for index,i in enumerate(location_of_alpha[2:]):
                ak=conformation_rotated[i,:]
                #print(ak)
                if not el.limitation(sign,a1_r_plane_norm_vec,distance,ak):
                    limit_ak=[False,(index+4)]
                    break
            conformation_rotated=np.matmul(conformation_rotated,rot_mat.T)
            #This line is for the debugging which will give out the residue that is forbiddened
            #result[current_frame,current_rotation]=limit_ak
            #This line will only give out the True or False value
            result[current_frame,current_rotation]=limit_ak[0]
            #print(limit_ak)
            #print(result[current_frame,current_rotation])
            current_rotation+=1
        current_frame+=1   
        #now = datetime.now().time()
        #print("now =", now)
        #print(result.shape)
    return result