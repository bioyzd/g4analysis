#*- coding : utf-8 -*-

'''
Created on 2011-01-18\n
The twist angle for G-DNA is defined from a JCTC article. which DOI is 10.1021/ct100253m.
The twist angle is defined using the angle between the line of C1' atoms in a G-quartet layer. 
@version: 0.1.0
@author: zhuh
@change:
    - 2011-01-18\n
        - Create this file.
    - 2011-01-19\n
        - finish the function B{Get_twist_in_GDNA()} and test it
    - 2011-01-25\n
        - modified the function B{Get_twist_in_GDNA()}, \
                so both gro and pdb can be used for coor_file.
        - add the version to B{0.1.0}
'''

import MDAnalysis
import numpy
import math
import PDB
import GRO
import Simple_atom
import usage

def Get_twist_in_GDNA(traj_file,coor_file,base_list_1,base_list_2,output_file):
    '''
    Input the layer 1 (G11,G12,G13,G14) and layer 2 (G21,G22,G23,G24),Calculate the angle 
    between G1i-G1(i+1) and G2i-G2(i+1). write the result to output_file.\n
    B{traj_file:} the GMX trajectory file, in trr or xtc format.\n
    B{coor_file:} the GMX coordinate file, in pdb or gro format.\n
    B{base_list_1:} the frist group contain four guanine bases.\n
    B{base_list_2:} the second group contain four guanine bases.\n
    B{output_file:} the output file. 
    '''
    C1_list_1=[]
    C1_list_2=[]
    print " init......"
    fp=open(output_file,"w")
    fp.write("#Group 1: ")
    for i in base_list_1:
        fp.write("%d\t " %i)
    fp.write("\n")
    fp.write("#Group 2: ")
    for i in base_list_2:
        fp.write("%d\t " %i)
    fp.write("\n")
    fp.write("#time\t angle_1 \t angle _2 \t angle _3 \t angle_4\n")
    Atom_list=[]
    if coor_file.endswith(".pdb"):
        Atom_list=PDB.Read_PDB_2_SimpleAtom_list(coor_file)
    elif coor_file.endswith(".gro"):
        Atom_list=GRO.Read_GRO_2_SimpleAtom_list(coor_file)
    else:
        pass
    for base in base_list_1:
        atom_list=Simple_atom.Get_Atom_in_residue(Atom_list,base)
        for atom in atom_list:
            if atom.atom_name =="C1'":
                C1_list_1.append(atom.atom_serial)

    for base in base_list_2:
        atom_list=Simple_atom.Get_Atom_in_residue(Atom_list,base)
        for atom in atom_list:
            if atom.atom_name =="C1'":
                C1_list_2.append(atom.atom_serial)

    u=MDAnalysis.Universe(coor_file,traj_file)    
    for ts in u.trajectory:
        angle=[]
        for i in range(4):
            vector1=[ts._x[C1_list_1[(i+1)%4]]-ts._x[C1_list_1[i]],\
                    ts._y[C1_list_1[(i+1)%4]]-ts._y[C1_list_1[i]],\
                    ts._z[C1_list_1[(i+1)%4]]-ts._z[C1_list_1[i]]]
            vector2=[ts._x[C1_list_2[(i+1)%4]]-ts._x[C1_list_2[i]],\
                    ts._y[C1_list_2[(i+1)%4]]-ts._y[C1_list_2[i]],\
                    ts._z[C1_list_2[(i+1)%4]]-ts._z[C1_list_2[i]]]
            vector1=numpy.array(vector1)
            vector2=numpy.array(vector2)
            gamma=numpy.dot(vector1,vector2)/(math.sqrt(numpy.dot(vector1,vector1)*numpy.dot(vector2,vector2)))
            angle.append(math.acos(gamma)/3.1416*180)

        fp.write("%6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f\n" \
                %(ts.time/1000,angle[0],angle[1],angle[2],angle[3]))             
        #if ts.frame % 100 ==0:
        #    print " analysis frame %6d......" %ts.frame
        usage.echo(" analysis frame %6d......\r" %ts.frame)
    fp.close()
    print "The result are in the file: ",output_file


def Get_twist_in_GDNA2(traj_file,coor_file,base_list_1,base_list_2,output_file):
    '''
    Input the layer 1 (G11,G12,G13,G14) and layer 2 (G21,G22,G23,G24),Calculate the angle 
    between G1i-G1(i+1) and G2i-G2(i+1). write the result to output_file.\n
    B{traj_file:} the GMX trajectory file, in trr or xtc format.\n
    B{coor_file:} the GMX coordinate file, in pdb or gro format.\n
    B{base_list_1:} the frist group contain four guanine bases.\n
    B{base_list_2:} the second group contain four guanine bases.\n
    B{output_file:} the output file. 
    '''
    C1_list_1=[]
    C1_list_2=[]
    print " init......"
    fp=open(output_file,"w")
    fp.write("#Group 1: ")
    for i in base_list_1:
        fp.write("%d\t " %i)
    fp.write("\n")
    fp.write("#Group 2: ")
    for i in base_list_2:
        fp.write("%d\t " %i)
    fp.write("\n")
    fp.write("#time\t angle_1 \t angle _2 \t angle _3 \t angle_4\n")
    Atom_list=[]
    if coor_file.endswith(".pdb"):
        Atom_list=PDB.Read_PDB_2_SimpleAtom_list(coor_file)
    elif coor_file.endswith(".gro"):
        Atom_list=GRO.Read_GRO_2_SimpleAtom_list(coor_file)
    else:
        pass
    for base in base_list_1:
        atom_list=Simple_atom.Get_Atom_in_residue(Atom_list,base)
        for atom in atom_list:
            if atom.atom_name =="C1'":
                C1_list_1.append(atom.atom_serial)

    for base in base_list_2:
        atom_list=Simple_atom.Get_Atom_in_residue(Atom_list,base)
        for atom in atom_list:
            if atom.atom_name =="C1'":
                C1_list_2.append(atom.atom_serial)

    u=MDAnalysis.Universe(coor_file,traj_file)    
    for ts in u.trajectory:
        angle=[]
        for i in range(4):
            vector1=[ts._x[C1_list_1[(i+1)%4]]-ts._x[C1_list_1[i]],\
                    ts._y[C1_list_1[(i+1)%4]]-ts._y[C1_list_1[i]],\
                    ts._z[C1_list_1[(i+1)%4]]-ts._z[C1_list_1[i]]]
            vector2=[ts._x[C1_list_2[(i+1)%4]]-ts._x[C1_list_2[i]],\
                    ts._y[C1_list_2[(i+1)%4]]-ts._y[C1_list_2[i]],\
                    ts._z[C1_list_2[(i+1)%4]]-ts._z[C1_list_2[i]]]
            vector1=numpy.array(vector1)
            vector2=numpy.array(vector2)
            gamma=numpy.dot(vector1,vector2)/(math.sqrt(numpy.dot(vector1,vector1)*numpy.dot(vector2,vector2)))
            angle.append(math.acos(gamma)/3.1416*180)

        fp.write("%6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f\n" \
                %(ts.time/1000,angle[0],angle[1],angle[2],angle[3]))             
        #if ts.frame % 100 ==0:
        #    print " analysis frame %6d......" %ts.frame
        usage.echo(" analysis frame %6d......\r" %ts.frame)
    fp.close()
    print "The result are in the file: ",output_file



