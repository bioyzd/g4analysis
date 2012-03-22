#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
Create on 2011-01-15.
用来测试各函数。
@change:\n
    - 2011.01.15
        - test the function B{Get_ratate_matrix()}
    - 2011.01.24
        - Add test_da.py to this file.

'''
import numpy
from G4Analysis import parallel_analysis
from G4Analysis import DNA_matrix
from G4Analysis import twist_in_GDNA
from G4Analysis import DNA_param
from G4Analysis import DNA_analysis
import unittest
from G4Analysis import Simple_atom

gro_file="confout.gro"
pdb_file="confout.pdb"
traj_file="trajout.xtc"
pqr_file="MOL.pqr"
idx_file="index.ndx"
amb_file="2GKU.top"


def test_Get_rotate_matrix():
    exper_coor=numpy.array([
        [ 11.417, -2.904, -4.880],
        [ 10.759, -1.995, -5.662],
        [ 11.469, -0.913, -5.867],
        [ 12.638, -1.108, -5.156],
        [ 13.759, -0.273, -5.036],
        [ 14.767, -0.848, -4.249],
        [ 14.663, -2.116, -3.719],
        [ 13.625, -2.934, -3.830],
        [ 12.625, -2.328, -4.545]
        ])
    base_name='G'
    A,B=parallel_analysis.Get_rotate_matrix(exper_coor,base_name)
    print A
    print B

def test_Get_RMSD():
    exper_coor=numpy.array([
        [ 11.417, -2.904, -4.880],
        [ 10.759, -1.995, -5.662],
        [ 11.469, -0.913, -5.867],
        [ 12.638, -1.108, -5.156],
        [ 13.759, -0.273, -5.036],
        [ 14.767, -0.848, -4.249],
        [ 14.663, -2.116, -3.719],
        [ 13.625, -2.934, -3.830],
        [ 12.625, -2.328, -4.545]
        ])
    base_name='G'
    rotation=numpy.array([
        [-0.2331, -0.8862, -0.4004],
        [ 0.8249, -0.3983,  0.4012],
        [-0.5150, -0.2368,  0.8238]
        ])
    origin=numpy.array([15.1632 , -0.0362 , -4.4678])
    A=parallel_analysis.Get_RMSD(exper_coor,base_name,origin,rotation)
    print A


def test_Rotate_2_vector():
    vector1=[-0.3781,0.5388,0.7528]
    vector2=[-0.4326,0.6579,0.6165]
    resu=DNA_matrix.Rotate_2_vector(vector1,vector2)
    print resu
    resu2=DNA_matrix.Rotate_2_vector(vector2,vector1)
    print resu2

def test_Get_parallel_result():
#    list1=[2]
#    list2=[3,9,17,21]
    list1=[[2,8,14,20]]
    list2=[[3,9,15,21]]
#    traj_file="/home/zhuh/Mount/99/ExPag/G4_ligand/2GKU-NMR/traj.xtc"
#    coor_file="/home/zhuh/Mount/99/ExPag/G4_ligand/2GKU-NMR/origin.pdb"
    global traj_file
    global gro_file
    parallel_analysis.Get_parallel_result(traj_file,gro_file,list1,list2,["parallel.xvg"])

def test_Get_parallel_fromTOP():
#    list1=[2]
#    list2=[3,9,17,21]
    list1=[2,8,14,20]
    list2=[3,9,15,21]
#    traj_file="/home/zhuh/Mount/99/ExPag/G4_ligand/2GKU-NMR/traj.xtc"
#    coor_file="/home/zhuh/Mount/99/ExPag/G4_ligand/2GKU-NMR/origin.pdb"
    global traj_file
    global gro_file
    parallel_analysis.Get_parallel_fromTOP(gro_file,list1,list2)



def test_Get_twist_in_GDNA():
    list1=[2,8,14,20]
    list2=[3,9,15,21]
    global traj_file
    global gro_file
    twist_in_GDNA.Get_twist_in_GDNA(traj_file,gro_file,list1,list2,"twist.xvg")

def test_Get_twist_in_GDNA2():
    list1=[[2,8,14,20]]
    list2=[[3,9,15,21]]
    global traj_file
    global gro_file
    twist_in_GDNA.Get_twist_in_GDNA2(traj_file,gro_file,list1,list2,["twist2.xvg"])

def test_middle_frame():
    A_G1_R=[[-0.2331,-0.8862,-0.4004],\
            [ 0.8249,-0.3983, 0.4012],\
            [-0.5150,-0.2368, 0.8238]]
    A_G1_O=[15.1632,-0.0362,-4.4678]
    B_C8_R=[[-0.2339,-0.9100,-0.3422],\
            [ 0.7496,-0.3930, 0.5326],\
            [-0.6191,-0.1320, 0.7741]]
    B_C8_O=[14.9142,0.2803,-4.7498]
    A_G2_R=[[-0.6807,-0.6274,-0.3781],\
            [ 0.3893,-0.7471, 0.5388],\
            [-0.6205, 0.2195, 0.7528]]
    A_G2_O=[14.8757,2.9250,-2.4635]
    B_C7_R=[[-0.5797,-0.6905,-0.4326],\
            [ 0.3207,-0.6814, 0.6579],\
            [-0.7491, 0.2426,0.6165]]
    B_C7_O=[14.4982,3.0313,-2.3001]
    middle_R_1,middle_O_1=DNA_param.middle_frame(A_G1_R,B_C8_R,A_G1_O,B_C8_O)
 #   print middle_R
 #   print middle_O
    middle_R_2,middle_O_2=DNA_param.middle_frame(A_G2_R,B_C7_R,A_G2_O,B_C7_O)
 #   print middle_R
 #   print middle_O
 #   print DNA_param.base_pair_parameters(A_G1_R,B_C8_R,A_G1_O,B_C8_O)
    print DNA_param.base_step_parameters(middle_R_1,middle_R_2,middle_O_1,middle_O_2)

def test_Get_pair_fromTOP():
    pdb_file="bdl084.pdb"

   
    print "-*"*35
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[1],[24],"pair")
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    print "-*"*35

    a=DNA_analysis.Get_para_fromTOP(pdb_file,[1,24],[2,23])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[2,23],[3,22])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[3,22],[4,21])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[4,21],[5,20])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[5,20],[6,19])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[6,19],[7,18])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[7,18],[8,17])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[8,17],[9,16])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[9,16],[10,15])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[10,15],[11,14])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   
    a=DNA_analysis.Get_para_fromTOP(pdb_file,[11,14],[12,13])
    print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])   



if __name__=="__main__":
#    test_Get_rotate_matrix()
#    test_Rotate_2_vector()
#    test_Get_parallel_result()
#    test_Get_parallel_fromTOP()
#    test_Get_RMSD()   
#    test_Get_twist_in_GDNA()
#    test_Get_twist_in_GDNA2()
#    test_GRO_Get_atom_List()
#    test_PDB_Get_atom_List()
#    test_Read_data_4_average()
#    test_Index_Atomlist_2_Index()
#    test_Traj_2_coor()
#    test_PQR_Get_Atom_list()
#    test_Index_Read_index_to_Inclass()
#    test_amber_top()
#    test_move_atom()
#    test_middle_frame()
    test_Get_pair_fromTOP()