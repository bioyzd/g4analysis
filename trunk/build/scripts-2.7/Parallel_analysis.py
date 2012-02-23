#!/usr/bin/python
'''
Create on 2011-01-19\n
@change:
    - 2011-01-20.\n
        - Now it's runing.
    - 2011-01-22.\n
        - Add code for checking the input. Make sure the user input 
        the right size and number of the base for a group.
    2011-02-01.\n
        - Modified this script, So both *.gro and *.pdb are allowd for \
                coordinare input file.
        - Changed the version to 0.2.1\n
    2011-04-26.\n 
        - Modified the usage part.
'''

import sys
import string
import re
import getopt
import os
from G4Analysis import Simple_atom
from G4Analysis import parallel_analysis
from G4Analysis import usage

def Usage(coor_file="coor_file",traj_file="traj_file",output_file="output_file",parm_file="para_analysis.in",skip=1,show_help="yes",calcu_rmsd="False"):
    '''
    print the usage information.
    '''
    print ""
    usage.File_input()
    usage.Print_line()
    usage.Coor_file(coor_file,"Input")
    usage.Traj_file(traj_file,"Input")
    usage.Xvgr_file(output_file,"Input")
    usage.Parm_file(parm_file,"Input")
    usage.Print_line()
    print ""
    usage.Type_input()
    usage.Print_line()
    usage.Show("--rmsd","bool",calcu_rmsd,"Calculate the RMSD of DNA bases groups.")
    usage.Show_skip(skip)
    usage.Show_help(show_help)
    print ""

def Check_argv(argv):
    '''
    Check the sys.argv list. and return a hash contain the input info.
    '''
    coor_file=""
    traj_file=""
    output_file=""
    parm_file=""
    skip=1
    resu_hash={}
    calcu_rmsd=False

    try:
        opts,args=getopt.getopt(sys.argv[1:],"p:f:o:i:h",["skip=","rmsd"])
    except getopt.GetoptError,e:
        print e
        sys.exit()

#    if len(opts) != 3 and len(opts) != 4 :
#        Usage()
#        sys.exit()

    for a,b in opts:
        if a == "-p" :
            coor_file=b
        elif a == "-f" :
            traj_file=b
        elif a == "-o":
            output_file=b
        elif a == "--skip" :
            try:
                skip = int(b)
            except:
                pass
        elif a == "--rmsd":
            calcu_rmsd=True
        elif a=="-i":
            parm_file=b
        elif a=="-h":
            Usage()
            sys.exit()


    if os.path.isfile(coor_file) and os.path.isfile(traj_file):
        resu_hash["coor_file"]=coor_file
        resu_hash["traj_file"]=traj_file
        resu_hash["parm_file"]=parm_file
        resu_hash["output_file"]=output_file

    if output_file=="" and parm_file=="":
        print "Error: No file for output."
        sys.exit()

    resu_hash["skip"]=skip
    resu_hash["calcu_rmsd"]=calcu_rmsd
    Usage(coor_file,traj_file,output_file,parm_file,skip,"no",calcu_rmsd)
    return resu_hash


def Print_ProgInfo():
    '''
    Print program information. the version, some notice and tips.
    '''
    print " "
    print "  "*12,"-)","Parallel_analysis","(-"
    print " "
    print "  "*12,"-)"," Version: %s " %usage.version ,"(-" 
    print " "
    Print_Description()
    print " Note:You should choose 1, 2 or 4 bases for a group."
    print " "

def Print_Description():
    '''
    Print the description of this script.
    '''
    print "DESCRIPTION"
    print "-----------"
    a='''Parallel_analysis.py is used for calculate the distance and angle between two bases groups.
usually a group contain 1, 2 or 4 bases in a plane.
The angle is useful to analysis the base stack. Two stack bases usually have a small angle
and fluctuation.
If the opition "--rmsd" used, only one bases group will be selected and the RMSD in z-axis 
for this group will be calculated.
Usage:
    interactive format:
    using -o result.xvg
    script format:
    using -i parm.in
para.in format:
    group_1(ID1:ID2:...) group_2(ID1:ID2:...) result.xvg
    '''
    print a
    print ""

if __name__=="__main__":
    Print_ProgInfo()
    argc=len(sys.argv)
    resu=Check_argv(sys.argv)
    is_get_dt=False
    have_parm_file=os.path.isfile(resu["parm_file"])

    list_group_1=list()
    list_group_2=list()
    list_output=list()


    if resu["traj_file"].endswith("mdcrd") or resu["traj_file"].endswith("dcd"):
        dt=raw_input("Input the time step between frames for Amber trajectory file (ps).")
        is_get_dt=True

    if resu["calcu_rmsd"]==False:
        if have_parm_file:
            fp=open(resu["parm_file"])
            lines=fp.readlines()
            for line in lines:
                try:
                    [group1,group2,outputname]=line.split()
                except ValueError,e:
                    print e
                    print "I guess you forget the --rmsd for calculating the rmsd."
                g1=group1.split(":")
                g2=group2.split(":")
                gg1=list()
                gg2=list()
                gg1=([int(i) for i in g1])
                gg2=([int(i) for i in g2])
                list_group_1.append(gg1)
                list_group_2.append(gg2)
                list_output.append(outputname)

        else:
            l1=Simple_atom.Get_list(resu["coor_file"],True)
            l2=Simple_atom.Get_list(resu["coor_file"],False)

            list_group_1.append(l1)
            list_group_2.append(l2)
            list_output.append(resu["output_file"])

        if is_get_dt:
            parallel_analysis.Get_parallel_result(resu["traj_file"],resu["coor_file"],list_group_1,list_group_2,list_output,resu["skip"],float(dt))
        else:
            parallel_analysis.Get_parallel_result(resu["traj_file"],resu["coor_file"],list_group_1,list_group_2,list_output,resu["skip"])

    else:
        if have_parm_file:
            fp=open(resu["parm_file"])
            lines=fp.readlines()
            for line in lines:
                [group1,outputname]=line.split()
                g1=group1.split(":")
                list_group_1.append([int(i) for i in g1])
                list_output.append(outputname)

        else:
            l1=Simple_atom.Get_list(resu["coor_file"],True)
            list_group_1.append(l1)

        if is_get_dt:
            parallel_analysis.Get_RMSD_result(resu["traj_file"],resu["coor_file"],list_group_1,list_output,resu["skip"],float(dt))
        else:
            parallel_analysis.Get_RMSD_result(resu["traj_file"],resu["coor_file"],list_group_1,list_output,resu["skip"])
