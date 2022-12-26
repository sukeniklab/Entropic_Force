# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 11:42:57 2022

Standard functions for entropic force analysis 
"""
# The input should be a np array containing the xyz coordinate

import numpy as np

# This function is to calculate the vector between 
# the alpha-carbons.
def vector(a1,a2):
    return (a2-a1)

# Magnitude of a vector
def maginitude(v):
    return np.linalg.norm(v)

# Inner product
def inner(v1,v2):
    return np.dot(v1,v2)

# Scalar projection 
# The first vector is our target and the second vector is the direction
def projection(v1,v2):
    return np.dot(v1,v2)/np.linalg.norm(v2)

# Angle between two vectors
def angel(v1,v2):
    return np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))

# When we write the actuall program, we can pass the vector norm to angel and 
# projection function to speed up our calculation
    
# Given a cone. The middle symmetry axis of the cone pass through a1,a2.
# The distance between apex and a1 is given by d. The angle betweenn the 
# generating line and the symmetry axis is also given by theta
# The d_a is the projection of other vector to the symmetry axis.
# d_a can be smaller than 1

def limit_distance(d, theta, d_a):
    if d_a<-d:
        return -999
    bot=(d+d_a)/np.tan(theta)
    return bot

# Following functions are for the entropic force rotation anaylsis
# Here we use a1a2 as a rotation axis and rotate the entire conformation to
# enhance sampling. 
# The direction and rotation will follow the Right-hand rule.

# Make a1 as the Coordinate system origin
# Minus a1 from every alpha carbon coordinate
def change_origin(mol,a1):
    return mol-a1

# Here we will find the normal vector of a plane to set a standard for
# rotation and find the constraint surface
def unit_normal_vector(v1,v2):
    normal_vec=np.cross(v1,v2)
    unit_normal_vec=normal_vec/np.linalg.norm(normal_vec)
    return unit_normal_vec
# The following rotation matrix describes the rotation around a certain 
# axis(vector) which goes through the origin. The rotation angle follows 
# the Right-hand rule. Source https://en.wikipedia.org/wiki/Rotation_matrix\
# The axis should be a unit vector. We will first normalize our vector and 
# then put them into this equation.

def rotation_matrix(angle,axis):
    axis=axis/np.linalg.norm(axis)
    cosangle=np.cos(angle)
    sinangle=np.sin(angle)
    return np.array([[cosangle+(1-cosangle)*np.square(axis[0]),(1-cosangle)*axis[0]*axis[1]-sinangle*axis[2],(1-cosangle)*axis[0]*axis[2]+sinangle*axis[1]],\
                            [(1-cosangle)*axis[1]*axis[0]+sinangle*axis[2],cosangle+(1-cosangle)*np.square(axis[1]),(1-cosangle)*axis[1]*axis[2]-sinangle*axis[0]],\
                            [(1-cosangle)*axis[2]*axis[0]-sinangle*axis[1],(1-cosangle)*axis[2]*axis[1]+sinangle*axis[0],cosangle+(1-cosangle)*np.square(axis[2])]])
# Here we want to find a vector a1a4. a1a4 is on the same plane as a1a2a3. 
# The angle between a1a2 and a1a4 is theta(given). There are two possible 
# a4 for each a1a2. To keep only one possible result, we will get a1a4 by
# rotating a1a2 by theta toward a1a3. This can be guranteed by using
# a1a2xa1a3 as axis and following the right hand rule.
# The axis should be a unit vector.
def vector_rotation(rotation_axis,vector,theta):
    matrix_r=rotation_matrix(theta,rotation_axis)
    vector_rotated=np.dot(matrix_r,vector)
    return vector_rotated    

# Based on the a1a4 and the normal vector of a1a2a3 (go through the a1), we 
# can determine a plane that is parallel to our target plane (go through the a1)
# a1_parallel_normal_vector = unit_normal_vector(a1a4,norm_vector)
# a1_paraller_normal_vector_unit= a1_parallel_normal_vector/np.linalg.norm(a1_parallel_normal_vector)
    
def distance_to_plane(unit_norm_vector,vector):
    return np.dot(unit_norm_vector,vector)

# Conformation rotation
def conformation_rotation(rotation_axis,conformation,theta):
    matrix_r=rotation_matrix(theta,rotation_axis)
    conformation_rotated=np.matmul(conformation,matrix_r.T)
    return conformation_rotated

# Determine whether the conformation is allowed
def limitation(plus_minus,normal_vector,distance,vector):
    v_distance=np.dot(vector,normal_vector)+plus_minus*distance
    #print(v_distance,vector,normal_vector,distance,plus_minus)
    if v_distance*plus_minus>0:
    # Allowed    
        return True
    else:
    # Prohibited
        return False

    