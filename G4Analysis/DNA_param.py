import numpy
import math
import DNA_matrix

def base_pair_parameters(rotation_matrix_1,rotation_matrix_2,origin_vector_1,origin_vector_2):

    y1_vector=[rotation_matrix_1[i][1] for i in range(3)]
    y2_vector=[rotation_matrix_2[i][1] for i in range(3)]

    gamma=math.acos(numpy.dot(y1_vector,y2_vector))  
    #calculate buckleopening angle.
    b0_vector=numpy.cross(y2_vector,y1_vector) 
    #buckle-opening axis
    b0_vector=b0_vector/math.sqrt(b0_vector.T * b0_vector)
    # normalize the b0 

    rotate_matrix_1=DNA_matrix.Rotate_matrix(-gamma/2,b0_vector)
    rotate_matrix_2=DNA_matrix.Rotate_matrix(gamma/2,b0_vector)
    #create rotate_matrix_2 

    transformed_orientation_matrix_1=numpy.matrix(rotate_matrix_1)*numpy.matirx(rotation_matrix_1)
    transformed_orientation_matrix_2=numpy.matrix(rotate_matrix_2)*numpy.matirx(rotation_matrix_2)
    #get the orientation matirx of transformed bases.

    MBT_matirx=(transformed_orientation_matrix_1+transformed_orientation_matrix_2)/2 
#//        MBT_origin_vector[i]=(origin_vector_1+origin_vector_2)/2 
   # // create the MBT matrix 
   
    for i in range(3):
        temp_vect=numpy.array([MBT_matirx[j,i] for j in range(3)])
        temp_vect=temp_vect/math.sqrt(temp_vect.T*temp_vect)
        for j in range(3):
            MBT_matirx[j,i]=temp_vect[j]
   #normalization the row.

    x_vector_temp_1=[transformed_orientation_matrix_1[i][0] for i in range(3)]
    x_vector_temp_2=[transformed_orientation_matrix_2[i][0] for i in range(3)]
    y_vector_MBT=[MBT_matirx[i][1] for i in range(3)]

    omega=math.acos(numpy.dot(x_vector_temp_1,x_vector_temp_2)) 
    cross_result_temp=numpy.cross(x_vector_temp_2,x_vector_temp_1) 

    if(numpy.dot(cross_result_temp,y_vector_MBT)<0):
        omega=-omega 

    x_vector_MBT=[MBT_matirx[i][0] for i in range(i)] 

    phi=math.acos(numpy.dot(b0_vector,x_vector_MBT)) 

    cross_b0_xmat_temp=numpy.cross(b0_vector,x_vector_MBT) 

    if(numpy.dot(cross_b0_xmat_temp,y_vector_MBT)<0):
        phi= - phi 

    kappa=gamma * math.cos(phi) 
    sigma=gamma * math.sin(phi) 

    #get kappa and sigma 

    displacement=numpy.array(origin_vector_1-origin_vector_2)*MBT_matirx

    shear=displacement[0] 
    stretch=displacement[1] 
    stagger=displacement[2] 
    buckle=kappa/PI * 180  
    propeller=omega /PI * 180 
    opening=sigma/PI * 180 

    return shear,stretch,stagger,buckle, propeller,opening


def base_step_paramaters(rotation_matrix_1,rotation_matrix_2,origin_vector_1,origin_vector_2):

    z1_vector=[rotation_matrix_1[i][2] for i in range(3)]
    z2_vector=[rotation_matrix_2[i][2] for i in range(3)]

    gamma=math.acos(numpy.dot(z1_vector,z2_vector,3)) 
    #   //calculate buckleopening angle.
    rt_vector=numpy.cross(z1_vector,z2_vector) 
    #//buckle-opening axis
    
    rt_vector=rt_vector/numpy.dot(rt_vector.T*rt_vector)
    #   //   normalize the rt_vector b0 

    rotate_matrix_1=DNA_matrix.Rotate_matrix(gamma/2,rt_vector)
    rotate_matrix_2=DNA_matrix.Rotate_matrix(-gamma/2,rt_vector)
    #//create rotate_matrix_1 
    #//create rotate_matrix_2 

    transformed_orientation_matrix_1=rotate_matrix_1*rotation_matrix_1 
    transformed_orientation_matrix_2=rotate_matrix_2*rotation_matrix_2 
    #//get the orientation matirx of transformed bases.

    MST_matirx=(transformed_orientation_matrix_1+transformed_orientation_matrix_2)/2 

    for i in range(3):
        temp_vect=numpy.array([MST_matirx[j,i] for j in range(3)])
        temp_vect=temp_vect/math.sqrt(temp_vect.T*temp_vect)
        for j in range(3):
            MST_matirx[j,i]=temp_vect[j]
   #normalization the row.

    y_vector_temp_1=[transformed_orientation_matrix_1[i][1] for i in range(i)]
    y_vector_temp_2=[transformed_orientation_matrix_2[i][1] for i in range(i)]
    z_vector_MST=[MST_matirx[i][2] for i in range(3)]

    omega=math.acos(numpy.dot(y_vector_temp_1,y_vector_temp_2)) 
    cross_result_temp=numpy.cross(y_vector_temp_1,y_vector_temp_2) 

    if(numpy.dot(cross_result_temp,z_vector_MST)<0):
        omega=-omega 
    # // for omega

    phi=math.acos(numpy.dot(rt_vector,y_vector_MST)) 

    cross_rt_ymat_temp=numpy.cross(rt_vector,y_vector_MST)

    if(numpy.dot(cross_rt_ymat_temp,z_vector_MST)<0):
        phi= - phi 
    #//for phi 

    kappa=gamma * math.cos(phi) 
    sigma=gamma * math.sin(phi) 

    #//get kappa and sigma 
    

    displacement=numpy.array(origin_vector_2[i] - origin_vector_1[i])*MST_matirx 

    shift=displacement[0] 
    slide=displacement[1] 
    rise=displacement[2] 
    roll=kappa/PI * 180  
    twist=omega /PI * 180 
    tilt=sigma/PI * 180 

    return shift,slide,rise,roll,twist,tilt


def middle_frame(rotation_matrix_1,rotation_matrix_2,origin_vector_1,origin_vector_2):

    y1_vector=[rotation_matrix_1[i][1] for i in range(3)]
    y2_vector=[rotation_matrix_2[i][1] for i in range(3)]

    gamma=math.acos(numpy.dot(y1_vector,y2_vector)) 
    #  //calculate buckleopening angle.
    b0_vector=numpy.cross(y2_vector,y1_vector)
    # //buckle-opening axis
#    print b0_vector
    
    b0_vector=b0_vector/math.sqrt(numpy.dot(b0_vector,b0_vector.T)) 
    #   //   normalize the b0 

    rotate_matrix_1=DNA_matrix.Rotate_matrix(-gamma/2,b0_vector)
    rotate_matrix_2=DNA_matrix.Rotate_matrix(gamma/2,b0_vector)
    #//create rotate_matrix_1 
    #//create rotate_matrix_2 

    transformed_orientation_matrix_1=numpy.matrix(rotate_matrix_1)*numpy.matrix(rotation_matrix_1) 
    transformed_orientation_matrix_1=DNA_matrix.Norm_matrix_in_row(transformed_orientation_matrix_1)
    
    transformed_orientation_matrix_2=numpy.matrix(rotate_matrix_2)*numpy.matrix(rotation_matrix_2) 
    transformed_orientation_matrix_2=DNA_matrix.Norm_matrix_in_row(transformed_orientation_matrix_2)
#//get the orientation matirx of transformed bases.

    middle_matirx=(transformed_orientation_matrix_1+transformed_orientation_matrix_2)/2 
#    print middle_matirx

    middle_matirx=DNA_matrix.Norm_matrix_in_row(middle_matirx)
   #normalization the row.

    middle_origin_vector=(numpy.array(origin_vector_1)+numpy.array(origin_vector_2))/2 

    return middle_matirx,middle_origin_vector

