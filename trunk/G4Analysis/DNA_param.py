import numpy
import math
import DNA_matrix

def base_pair_parameters(rotation_1,rotation_2,origin_1,origin_2):

    y1_vector=numpy.array([rotation_1[i][1] for i in range(3)])
    y2_vector=numpy.array([rotation_2[i][1] for i in range(3)])

    gamma=math.acos(numpy.dot(y1_vector,y2_vector))  
    #calculate buckleopening angle.
    b0_vector=numpy.cross(y2_vector,y1_vector) 
    #buckle-opening axis
    b0_vector=b0_vector/math.sqrt(numpy.dot(b0_vector,b0_vector.T))
    # normalize the b0 

    rotate_matrix_1=DNA_matrix.Rotate_matrix(-gamma/2,b0_vector)
    rotate_matrix_2=DNA_matrix.Rotate_matrix(gamma/2,b0_vector)
    #create rotate_matrix_2 

    trans_orient_1=numpy.matrix(rotate_matrix_1)*numpy.matrix(rotation_1)
    trans_orient_2=numpy.matrix(rotate_matrix_2)*numpy.matrix(rotation_2)
    #get the orientation matirx of transformed bases.

    MBT_matirx=(trans_orient_1+trans_orient_2)/2 
#//        MBT_origin[i]=(origin_1+origin_2)/2 
   
    MBT_matirx=DNA_matrix.Norm_matrix_in_row(MBT_matirx)
   #normalization the row.

    x_vector_temp_1=numpy.array([trans_orient_1[i,0] for i in range(3)])
    x_vector_temp_2=numpy.array([trans_orient_2[i,0] for i in range(3)])
#    print x_vector_temp_1
#    print x_vector_temp_2
    x_vector_temp_1=x_vector_temp_1/math.sqrt(numpy.dot(x_vector_temp_1,x_vector_temp_1.T))
    x_vector_temp_2=x_vector_temp_2/math.sqrt(numpy.dot(x_vector_temp_2,x_vector_temp_2.T))
#    print x_vector_temp_1
#    print x_vector_temp_2
    y_vector_MBT=numpy.array([MBT_matirx[i,1] for i in range(3)])

    omega=math.acos(numpy.dot(x_vector_temp_1,x_vector_temp_2.T)) 
    cross_result_temp=numpy.cross(x_vector_temp_2,x_vector_temp_1) 

    if(numpy.dot(cross_result_temp,y_vector_MBT.T)<0):
        omega=-omega 

    x_vector_MBT=numpy.array([MBT_matirx[i,0] for i in range(3)])

    phi=math.acos(numpy.dot(b0_vector,x_vector_MBT.T)) 

    cross_b0_xmat_temp=numpy.cross(b0_vector,x_vector_MBT) 

    if(numpy.dot(cross_b0_xmat_temp,y_vector_MBT)<0):
        phi= - phi 

    kappa=gamma * math.cos(phi) 
    sigma=gamma * math.sin(phi) 

    #get kappa and sigma 

    displacement=numpy.matrix(numpy.array(origin_1)-numpy.array(origin_2))*MBT_matirx

    shear=displacement[0,0] 
    stretch=displacement[0,1] 
    stagger=displacement[0,2] 
    buckle=kappa/math.pi * 180  
    propeller=omega /math.pi * 180 
    opening=sigma/math.pi * 180 

    return shear,stretch,stagger,buckle, propeller,opening

def base_step_parameters(rotation_1,rotation_2,origin_1,origin_2):

    z1_vector=numpy.array([rotation_1[i][2] for i in range(3)])
    z2_vector=numpy.array([rotation_2[i][2] for i in range(3)])

    gamma=math.acos(numpy.dot(z1_vector,z2_vector.T)) 
#    print gamma
    #   //calculate buckleopening angle.
    rt_vector=numpy.cross(z1_vector,z2_vector) 
    #//buckle-opening axis
    
    rt_vector=rt_vector/math.sqrt(numpy.dot(rt_vector,rt_vector.T))
    #   //   normalize the rt_vector b0 

    rotate_matrix_1=DNA_matrix.Rotate_matrix(gamma/2,rt_vector)
    rotate_matrix_2=DNA_matrix.Rotate_matrix(-gamma/2,rt_vector)
    #//create rotate_matrix_1, create rotate_matrix_2 

    trans_orient_1=numpy.matrix(rotate_matrix_1)*numpy.matrix(rotation_1)
    trans_orient_2=numpy.matrix(rotate_matrix_2)*numpy.matrix(rotation_2) 
#    print trans_orient_1
#    print trans_orient_2
    #//get the orientation matirx of transformed bases.

    MST_matirx=(trans_orient_1+trans_orient_2)/2 
    MST_matirx=DNA_matrix.Norm_matrix_in_row(MST_matirx)
    #normalization the row.

    y_vector_temp_1=numpy.array([trans_orient_1[i,1] for i in range(3)])
    y_vector_temp_2=numpy.array([trans_orient_2[i,1] for i in range(3)])
    y_vector_temp_1=y_vector_temp_1/math.sqrt(numpy.dot(y_vector_temp_1,y_vector_temp_1.T))
    y_vector_temp_2=y_vector_temp_2/math.sqrt(numpy.dot(y_vector_temp_2,y_vector_temp_2.T))
    z_vector_MST   =numpy.array([MST_matirx[i,2] for i in range(3)])
    y_vector_MST   =numpy.array([MST_matirx[i,1] for i in range(3)])

    omega=math.acos(numpy.dot(y_vector_temp_1,y_vector_temp_2.T)) 
    cross_result_temp=numpy.cross(y_vector_temp_1,y_vector_temp_2) 

    if(numpy.dot(cross_result_temp,z_vector_MST)<0):
        omega=-omega 
    # // for omega

    phi=math.acos(numpy.dot(rt_vector,y_vector_MST.T)) 

    cross_rt_ymat_temp=numpy.cross(rt_vector,y_vector_MST)

    if(numpy.dot(cross_rt_ymat_temp,z_vector_MST)<0):
        phi= - phi 
    #//for phi 

    kappa=gamma * math.cos(phi) 
    sigma=gamma * math.sin(phi) 

    #//get kappa and sigma 

    displacement=numpy.matrix(origin_2 - origin_1)*MST_matirx 

    shift=displacement[0,0] 
    slide=displacement[0,1] 
    rise=displacement[0,2] 
    roll=kappa/math.pi * 180  
    twist=omega /math.pi * 180 
    tilt=sigma/math.pi * 180 

    return shift,slide,rise,roll,twist,tilt


def middle_frame(rotation_1,rotation_2,origin_1,origin_2):

    y1_vector=[rotation_1[i][1] for i in range(3)]
    y2_vector=[rotation_2[i][1] for i in range(3)]

    gamma=math.acos(numpy.dot(y1_vector,y2_vector)) 
    #  //calculate buckleopening angle.
    b0_vector=numpy.cross(y2_vector,y1_vector)
    # //buckle-opening axis
    
    b0_vector=b0_vector/math.sqrt(numpy.dot(b0_vector,b0_vector.T)) 
    #   //   normalize the b0 

    rotate_matrix_1=DNA_matrix.Rotate_matrix(-gamma/2,b0_vector)
    rotate_matrix_2=DNA_matrix.Rotate_matrix(gamma/2,b0_vector)
    #//create rotate_matrix_1 
    #//create rotate_matrix_2 

    trans_orient_1=numpy.matrix(rotate_matrix_1)*numpy.matrix(rotation_1) 
    trans_orient_1=DNA_matrix.Norm_matrix_in_row(trans_orient_1)
    
    trans_orient_2=numpy.matrix(rotate_matrix_2)*numpy.matrix(rotation_2) 
    trans_orient_2=DNA_matrix.Norm_matrix_in_row(trans_orient_2)
#//get the orientation matirx of transformed bases.

    middle_matirx=(trans_orient_1+trans_orient_2)/2 

    middle_matirx=DNA_matrix.Norm_matrix_in_row(middle_matirx)
   #normalization the row.

    middle_origin=(numpy.array(origin_1)+numpy.array(origin_2))/2 

    return numpy.array(middle_matirx),middle_origin

