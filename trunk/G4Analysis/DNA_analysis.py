#1/usr/bin/env python
#-*- coding:utf-8 -*-
'''
Created at 2012.03.06
'''

import Simple_atom
import time as Time
import sys
import os
import usage
import MDAnalysis
import DNA_param
from math import cos, sin, sqrt
import math
import numpy
from numpy import matrix
from numpy import dot

import DNA_matrix

def Get_parm_fromTRJ(traj_file, coor_file, base_list_1, base_list_2, output_name,CALCU="step",skip=1, dt=1,begin=0,end=-1):
    '''
    Note: this function not finish yet.
    base_list_1=[[base or base pair],[base or base pair],...]
    base_list_2=[[base or base pair],[base or base pair],...]
    '''
    print "  init......"
#    if coor_file.endswith('.top'):
#        Atom_list=amber_top
    START_TIME=Time.time()
    if len(base_list_1)==len(base_list_2):
        LIST_NUM=len(base_list_1)
    else:
        print "ERROR: The length of the base list not match."
        return -1
    if len(output_name)!=LIST_NUM:
        print "ERROR: The number of the output file not match the size of the base list."

    Atom_list=Simple_atom.Get_Simple_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)

    base_name_list_1=list()
    base_name_list_2=list()
    base_atom_list_1=list()
    base_atom_list_2=list()

    for i in range(LIST_NUM):

        if os.path.isfile(output_name[i]):
            print "backup %s to %s" %(output_name[i],"#"+output_name[i]+"#")
            try:
                os.rename(output_name[i],"#"+output_name[i]+"#")
            except OSError,e: 
                print e
                print "the file %s will be overwrited!" %output_name[i]

        fp = open(output_name[i], 'w')
        fp.write("#Group 1: ")
        for j in base_list_1[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#Group 2: ")
        for j in base_list_2[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#skip:%d\n" %skip)
        if CALCU=="step":
            fp.write("#  shift       slide        rise        tilt        roll       twist\n")
        else:
            fp.write("#  shear     stretch     stagger      buckle   propeller     opening\n")
        fp.close()

#        base_name_list_1.append( [residue_list[j-1] for j in base_list_1[i]])
#        base_name_list_2.append( [residue_list[j-1] for j in base_list_2[i]])

        temp_list=list()
        for m in base_list_1[i]:
            for n in residue_list:
                if n[1]==m:
                    temp_list.append(n)
                    break
                else:
                    pass
        base_name_list_1.append(temp_list)

        temp_list=list()
        for m in base_list_2[i]:
            for n in residue_list:
                if n[1]==m:
                    temp_list.append(n)
                    break
                else:
                    pass
        base_name_list_2.append(temp_list)


        base_atom_list_1.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_1[i]])
        base_atom_list_2.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_2[i]])


    u=MDAnalysis.Universe(coor_file,traj_file)

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        try:
            dt=u.trajectory.dt
        except:
            dt=0.0

    for ts in u.trajectory:
        time=float((ts.frame-1)*dt)
        if dt > 0.0:
            if time < float(begin):
                continue
            if time > float(end) and end !=-1:
                break

        if ts.frame % skip == 0 :
            for i in range(LIST_NUM):
                r1=[]
                '''the group 1 rotate list'''
                r2=[]
                '''the group 2 rotate list'''
                c1=[]
                '''the group 1 coordinate list'''
                c2=[]
                '''the group 2 coordinate list'''
                for m in range(len(base_name_list_1[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_1[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_1[i][m][0])
                    #base_name_list_1[index of the groups][index of the base of group 1][base_name,base_serial]
                    r1.append(result[0])
                    c1.append(result[1])

                for m in range(len(base_name_list_2[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_2[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[i][m][0])
                    r2.append(result[0])
                    c2.append(result[1])

                fp = open(output_name[i], 'a')

                if CALCU=="pair":
                    a=DNA_param.base_pair_parameters(r1[0],r2[0],c1[0],c2[0])
                    fp.write("%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n" %(a[0],a[1],a[2],a[3],a[4],a[5]))
                else:
                    middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
                    middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])
                    a=DNA_param.base_step_parameters(middle_r1,middle_r2,middle_c1,middle_c2)
                    fp.write("%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n" %(a[0],a[1],a[2],a[3],a[4],a[5]))

                if ts.frame % 10 ==0 and i==0:
                    NOW_TIME=Time.time()
                    usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))

                fp.close()

    print "The DNA helical analysis finished"
    print "The result are in file: %s" %output_name

def Get_para_fromTOP( coor_file, base_list_1, base_list_2,CALCU="step",PRINT=True):
    '''
    get pair parameters from coor_file
    '''

    Atom_list=Simple_atom.Get_Simple_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)
    atom_list=Simple_atom.Get_atom_list(coor_file)

    base_name_list_1=list()
    base_name_list_2=list()
    base_atom_list_1=list()
    base_atom_list_2=list()

    for m in base_list_1:
        for n in residue_list:
            if n[1]==m:
                base_name_list_1.append(n)
                break
            else:
                pass

    for m in base_list_2:
        for n in residue_list:
            if n[1]==m:
                base_name_list_2.append(n)
                break
            else:
                pass

    base_atom_list_1=[DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_1]
    base_atom_list_2=[DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_2]

    r1=[]
    '''the group 1 rotate list'''
    r2=[]
    '''the group 2 rotate list'''
    c1=[]
    '''the group 1 coordinate list'''
    c2=[]
    '''the group 2 coordinate list'''
    for m in range(len(base_name_list_1)):
        temp_list = [ [atom_list[x].atom_coor_x*10, atom_list[x].atom_coor_y*10,atom_list[x].atom_coor_z*10] \
                for x in base_atom_list_1[m] ]
#        print temp_list
        result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_1[m][0])
        r1.append(result[0])
        c1.append(result[1])

    for m in range(len(base_name_list_2)):
        temp_list = [ [atom_list[x].atom_coor_x*10, atom_list[x].atom_coor_y*10,atom_list[x].atom_coor_z*10] \
                for x in base_atom_list_2[m] ]
#        print temp_list
        result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[m][0])
        r2.append(result[0])
        c2.append(result[1])

    if CALCU=="pair":
        a=DNA_param.base_pair_parameters(r1[0],r2[0],c1[0],c2[0])
        if PRINT:
            print "   shear     stretch     stagger      buckle   propeller     opening"
            print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])
        return a
    else:
        middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
        middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])
        a=DNA_param.base_step_parameters(middle_r1,middle_r2,middle_c1,middle_c2)
        if PRINT:
            print "   shift       slide        rise        tilt        roll       twist"
            print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])
        return a


def Get_Dihedral_fromTOP( coor_file, base_list_1,PRINT=False):
    '''
    calculate the dihedral from the top file.

    '''
    Atom_list=Simple_atom.Get_atom_list(coor_file)

    for m in base_list_1:
        for i in Atom_list:
      #      print Atom_list[i].atom_name
            if Atom_list[i].atom_name in ["O3'","O3*"] and Atom_list[i].residue_serial==m-1: 
                index_O3_0=i

            elif Atom_list[i].atom_name=="P"   and Atom_list[i].residue_serial==m:
                index_P=i

            elif Atom_list[i].atom_name in ["O5'","O5*"] and Atom_list[i].residue_serial==m:
                index_O5=i
            elif Atom_list[i].atom_name in ["C5'","C5*"] and Atom_list[i].residue_serial==m:
                index_C5=i
            elif Atom_list[i].atom_name in ["C4'","C4*"] and Atom_list[i].residue_serial==m:
                index_C4=i
            elif Atom_list[i].atom_name in ["C3'","C3*"] and Atom_list[i].residue_serial==m:
                index_C3=i
            elif Atom_list[i].atom_name in ["O3'","O3*"] and Atom_list[i].residue_serial==m:
                index_O3=i

            elif Atom_list[i].atom_name in ["P"] and Atom_list[i].residue_serial==m+1:
                index_P_2=i
            elif Atom_list[i].atom_name  in ["O5'","O5*"] and Atom_list[i].residue_serial==m+1:
                index_O5_2=i

            elif Atom_list[i].atom_name in ["O4'","O4*"] and Atom_list[i].residue_serial==m:
                index_O4=i

            elif Atom_list[i].atom_name in ["C1'","C1*"] and Atom_list[i].residue_serial==m:
                index_C1=i

            elif Atom_list[i].atom_name in ["N9"] and Atom_list[i].residue_serial==m and \
                    ("A" in Atom_list[i].residue_name or "G" in Atom_list[i].residue_name):
                index_N_base=i
            elif Atom_list[i].atom_name=="N1" and Atom_list[i].residue_serial==m and \
                    ("C" in Atom_list[i].residue_name or "T" in Atom_list[i].residue_name):
                index_N_base=i
            elif Atom_list[i].atom_name=="C4" and Atom_list[i].residue_serial==m and \
                    ("A" in Atom_list[i].residue_name or "G" in Atom_list[i].residue_name):
                index_C_base=i
            elif Atom_list[i].atom_name=="C2" and Atom_list[i].residue_serial==m and \
                    ("C" in Atom_list[i].residue_name or "T" in Atom_list[i].residue_name):
                index_C_base=i

        alpha =DNA_param.Get_dihedral(Atom_list[index_O3_0],Atom_list[index_P ],Atom_list[index_O5 ],Atom_list[index_C5  ]) *180/math.pi
        beta  =DNA_param.Get_dihedral(Atom_list[index_P   ],Atom_list[index_O5],Atom_list[index_C5 ],Atom_list[index_C4  ]) *180/math.pi
        gamma =DNA_param.Get_dihedral(Atom_list[index_O5  ],Atom_list[index_C5],Atom_list[index_C4 ],Atom_list[index_C3  ]) *180/math.pi
        delta =DNA_param.Get_dihedral(Atom_list[index_C5  ],Atom_list[index_C4],Atom_list[index_C3 ],Atom_list[index_O3  ]) *180/math.pi
        epslon=DNA_param.Get_dihedral(Atom_list[index_C4  ],Atom_list[index_C3],Atom_list[index_O3 ],Atom_list[index_P_2 ]) *180/math.pi
        zeta  =DNA_param.Get_dihedral(Atom_list[index_C3  ],Atom_list[index_O3],Atom_list[index_P_2],Atom_list[index_O5_2]) *180/math.pi

        chi   =DNA_param.Get_dihedral(Atom_list[index_O4],Atom_list[index_C1],Atom_list[index_N_base],Atom_list[index_C_base]) *180/math.pi
        if PRINT == True:
            print "%10s%10s%10s%10s%10s%10s" %("alpha","beta","gamma","delta","epslon","zeta")
            print "%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f" %(alpha,beta,gamma,delta,epslon,zeta)

        return [alpha,beta,gamma,delta,epslon,zeta,chi]


Get_Dihedral_fromTRJ(traj_file, coor_file, base_list, output_name,skip=1, dt=1,begin=0,end=-1):
    START_TIME=Time.time()
    if len(base_list_1)==len(base_list_2):
        LIST_NUM=len(base_list_1)
    else:
        print "ERROR: The length of the base list not match."
        return -1
    if len(output_name)!=LIST_NUM:
        print "ERROR: The number of the output file not match the size of the base list."

    Atom_list=Simple_atom.Get_Simple_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)

    base_name_list_1=list()
    base_name_list_2=list()
    base_atom_list_1=list()
    base_atom_list_2=list()

    for i in range(LIST_NUM):

        if os.path.isfile(output_name[i]):
            print "backup %s to %s" %(output_name[i],"#"+output_name[i]+"#")
            try:
                os.rename(output_name[i],"#"+output_name[i]+"#")
            except OSError,e: 
                print e
                print "the file %s will be overwrited!" %output_name[i]

        fp = open(output_name[i], 'w')
        fp.write("#Group 1: ")
        for j in base_list_1[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#Group 2: ")
        for j in base_list_2[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#skip:%d\n" %skip)
        if CALCU=="step":
            fp.write("#  shift       slide        rise        tilt        roll       twist\n")
        else:
            fp.write("#  shear     stretch     stagger      buckle   propeller     opening\n")
        fp.close()

#        base_name_list_1.append( [residue_list[j-1] for j in base_list_1[i]])
#        base_name_list_2.append( [residue_list[j-1] for j in base_list_2[i]])

        temp_list=list()
        for m in base_list_1[i]:
            for n in residue_list:
                if n[1]==m:
                    temp_list.append(n)
                    break
                else:
                    pass
        base_name_list_1.append(temp_list)

        temp_list=list()
        for m in base_list_2[i]:
            for n in residue_list:
                if n[1]==m:
                    temp_list.append(n)
                    break
                else:
                    pass
        base_name_list_2.append(temp_list)


        base_atom_list_1.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_1[i]])
        base_atom_list_2.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_2[i]])


    u=MDAnalysis.Universe(coor_file,traj_file)

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        try:
            dt=u.trajectory.dt
        except:
            dt=0.0

    for ts in u.trajectory:
        time=float((ts.frame-1)*dt)
        if dt > 0.0:
            if time < float(begin):
                continue
            if time > float(end) and end !=-1:
                break

        if ts.frame % skip == 0 :
            for i in range(LIST_NUM):
                r1=[]
                '''the group 1 rotate list'''
                r2=[]
                '''the group 2 rotate list'''
                c1=[]
                '''the group 1 coordinate list'''
                c2=[]
                '''the group 2 coordinate list'''
                for m in range(len(base_name_list_1[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_1[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_1[i][m][0])
                    #base_name_list_1[index of the groups][index of the base of group 1][base_name,base_serial]
                    r1.append(result[0])
                    c1.append(result[1])

                for m in range(len(base_name_list_2[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_2[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[i][m][0])
                    r2.append(result[0])
                    c2.append(result[1])

                fp = open(output_name[i], 'a')

                if CALCU=="pair":
                    a=DNA_param.base_pair_parameters(r1[0],r2[0],c1[0],c2[0])
                    fp.write("%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n" %(a[0],a[1],a[2],a[3],a[4],a[5]))
                else:
                    middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
                    middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])
                    a=DNA_param.base_step_parameters(middle_r1,middle_r2,middle_c1,middle_c2)
                    fp.write("%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n" %(a[0],a[1],a[2],a[3],a[4],a[5]))

                if ts.frame % 10 ==0 and i==0:
                    NOW_TIME=Time.time()
                    usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))

                fp.close()

    print "The DNA helical analysis finished"
    print "The result are in file: %s" %output_name


