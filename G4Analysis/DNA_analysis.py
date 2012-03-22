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
    Reading the traj file and the coordinate file like *.pdb or *.gro. With the base serial choosed,  get the
    rotation matrix for this base. and write it's to a output file with some syntax.
    base_list_1 format list(list())
    base_list_2 format list(list())
    output_name format list()
    '''
    print "  init......"
#    if coor_file.endswith('.top'):
#        Atom_list=amber_top
    START_TIME=Time.time()
    LIST_NUM=len(base_list_1)

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
        fp.write("#time(ns)   distance(A)\t   angle(degree)    RMSD1(A)\t    RMSD2(A)\n")
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
                    c1.append(numpy.array(temp_list))
                    r1.append(result)

                for m in range(len(base_name_list_2[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_2[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[i][m][0])
                    c2.append(numpy.array(temp_list))
                    r2.append(result)

                orient_group_1,origin_group_1 = DNA_matrix.Get_group_rotmat(r1,len(base_name_list_1[i]))
                orient_group_2,origin_group_2 = DNA_matrix.Get_group_rotmat(r2,len(base_name_list_2[i]))
                RMSD1=DNA_matrix.Get_group_RMSD(base_name_list_1[i],c1,origin_group_1,orient_group_1)
                RMSD2=DNA_matrix.Get_group_RMSD(base_name_list_2[i],c2,origin_group_2,orient_group_2)
    
                if numpy.dot(orient_group_1, orient_group_2)<0:
                    orient_group_2 = orient_group_2*(-1)

                orient_total = DNA_matrix.Rotate_2_vector(orient_group_1, orient_group_2)
                dist_vector = numpy.array(origin_group_2)-numpy.array(origin_group_1)

                dist = abs(dot(orient_total, dist_vector))
                gamma = dot(orient_group_1, orient_group_2) 
                gamma = math.acos(gamma)*180/3.1416


                if ts.frame % 10 ==0 and i==0:
                    NOW_TIME=Time.time()
                    usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))

                fp = open(output_name[i], 'a')
                fp.write( " %7.4f\t  %6.3f\t   %6.3f\t    %6.4f\t   %6.4f\n" %(time/1000, dist, gamma,RMSD1,RMSD2))
                fp.close()

    print "The parallel analysis finished"
    print "The result are in file: %s" %output_name

def Get_para_fromTOP( coor_file, base_list_1, base_list_2,CALCU="step"):
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
        return DNA_param.base_pair_parameters(r1[0],r2[0],c1[0],c2[0])
    else:
        middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
        middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])
        return  DNA_param.base_step_parameters(middle_r1,middle_r2,middle_c1,middle_c2)

