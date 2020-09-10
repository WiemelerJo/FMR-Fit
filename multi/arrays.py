parameter_Table_names_L = ['dB','R','A']
parameter_table_names_L_final = []
parameter_Table_names_D = ['alpha','dB','R','A']
parameter_table_names_D_final = []

#Arrays to set the default values in SpinBox. #For linear: [slope, offset]
                                             #For L: [dB1, R1, A1, dB2, R2, A2, ...]
                                             #For D: [alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, ....]
default_linear = [0.5,1]
default_values_L = [0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7, 0.03,0.08,7]
default_values_D = [0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7, 0.0001,0.03,0.08,7]

#Arrays for default Boundaries. For Linear : [slope_min,slope_max, offset_min,offset_max] 
                               #For L: [dB1_min,dB1_max, R1_min,R1_max, ...]
                               #For D: [alpha1_min,alpha1_max, dB1_min,dB1_max ...]
default_boundaries_linear_min = [-5, -5]
default_boundaries_linear_max = [5, 5]
default_boundaries_L_min = [0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0, 0.0001,0,0]
default_boundaries_L_max = [0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15,0.5,0.15,15]
default_boundaries_D_min = [0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001, 0.0001,0.0001,0.0001,0.0001]
default_boundaries_D_max = [1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15, 1,0.5,0.15,15]

#Arrays for step size for spinbox Linear: [slope, offset]
                                 #L: [dB, R, A]
                                 #D: [alpha,dB, R, A]
default_stepsize_linear = [0.1, 0.1]
default_stepsize_L = [0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1, 0.0001,0.01,1]
default_stepsize_D = [0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1, 0.01,0.0001,0.01,1]

#Arrays for maximum boundary for spinbox Linear: [slope, offset]
                                             #L: [dB, R, A]
                                             #D: [alpha,dB, R, A]
default_maximum_bound_spinbox_linear = [1000,1000]
default_maximum_bound_spinbox_L = [1.5,3.0,15000, 1.5,3.0,15000, 1.5,3.0,15000, 1.5,3.0,15000, 1.5,3.0,15000, 1.5,3.0,15000, 1.5,3.0,15000, 1.5,3.0,15000, 1.5,3.0,15000, 1.5,3.0,15000]
default_maximum_bound_spinbox_D = [1.0,1.5,3.0,15000 ,1.0,1.5,3.0,15000 ,1.0,1.5,3.0,15000 ,1.0,1.5,3.0,15000 ,1.0,1.5,3.0,15000 ,1.0,1.5,3.0,15000 ,1.0,1.5,3.0,15000 ,1.0,1.5,3.0,15000 ,1.0,1.5,3.0,15000 ,1.0,1.5,3.0,15000]
