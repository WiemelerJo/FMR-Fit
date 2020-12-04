import numpy as np
# Gen arrays to set the default parameters in SpinBox.

# Arrays for names of parameters. Will be displayed at the side of the table
parameter_Table_names_L = ['dB', 'R', 'A']
parameter_Table_names_D = ['alpha', 'dB', 'R', 'A']

class Gen_array():
    def __init__(self,fit_num,**kws):
        self.increment = 0.001
        for key, val in kws.items():
            setattr(self, key, val)
        if not kws:
            self.increase = False
        self.increment_orig = self.increment
        self.fit_num = fit_num

        self.do_Stuff()

    def gen_array(self,ary:list,increase=False,increment=0.001):
        # ary now can now be any array with size for one function
        array = []
        if increase:
            for num in range(1,self.fit_num + 1):
                for i in ary:
                    array.append(i + increment)
                    increment += self.increment_orig
        else:
            for num in range(1,self.fit_num + 1):
                array.append(ary)
        array = np.asarray(array).flatten()
        return array

    def do_Stuff(self):
        # ---------------Create Arrays for the values according to the order in that they are created-------------------

        # Arrays to set the initial VALUES in SpinBox. #For linear: [slope, offset]
        # For L: [dB1, R1, A1, dB2, R2, A2, ...], For one func [0.03,0.08,7]
        # For D: [alpha1, dB1, R1, A1, alpha2, dB2, R2, A2, ....]  For one func [0.0001,0.03,0.08,7]
        self.default_linear = [0.5, 1]
        self.default_values_L = self.gen_array([0.03, 0.08, 7],self.increase,self.increment)
        self.default_values_D = self.gen_array([0.0001, 0.03, 0.08, 7],self.increase,self.increment)

        # Arrays for default Boundaries. For Linear : [slope_min,slope_max, offset_min,offset_max]
        # For L_min : [dB1_min,R1_min,A1_min]
        # For L_max : [dB1_max,R1_max,A1_max]
        # For D_min: [alpha1_min, dB1_min,...]
        # For D_max: [alpha1_max, dB1_max ...]
        self.default_boundaries_linear_min = [-5, -5]
        self.default_boundaries_linear_max = [5, 5]
        self.default_boundaries_L_min = self.gen_array([0.0001, 0, 0])
        self.default_boundaries_L_max = self.gen_array([0.5, 0.15, 15])
        self.default_boundaries_D_min = self.gen_array([0.0001, 0.0001, 0.0001, 0.0001])
        self.default_boundaries_D_max = self.gen_array([1, 0.5, 0.15, 15])

        # Arrays for step size for spinbox Linear: [slope, offset]
        # L: [dB, R, A]
        # D: [alpha,dB, R, A]
        self.default_stepsize_linear = [0.1, 0.1]
        self.default_stepsize_L = self.gen_array([0.0001, 0.01, 1])
        self.default_stepsize_D = self.gen_array([0.01, 0.0001, 0.01, 1])

        # Arrays for maximum boundary for spinbox Linear: [slope, offset]
        # L: [dB, R, A]
        # D: [alpha,dB, R, A]
        self.default_maximum_bound_spinbox_linear = [1000, 1000]
        self.default_maximum_bound_spinbox_L = self.gen_array([1.5, 3.0, 15000])
        self.default_maximum_bound_spinbox_D = self.gen_array([1.0, 1.5, 3.0, 15000])

#--------Usage---------------
#obj = Gen_array(2) # Define obj, with number of fit functions (fit_num)
#print(obj.default_maximum_bound_spinbox_D) # Adress Array of interest

#test = Gen_array(2,increase=True)