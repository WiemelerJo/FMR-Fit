import matplotlib as mp
mp.use('QT5Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
import math as m

class Params_Plot(object):
    """docstring for Params_Plot"""
    def __init__(self, params_fname,BOOL,index_model):
        #QThread.__init__(self)
        super().__init__()
        #self.params_fname = params_fname
        if BOOL == 'eigen':
            self.plot_eigen(params_fname,index_model)
        else:
            self.plot(params_fname,index_model)

    def plot(self,params_fname,index_model):

        #Todo:  Try to implement a method where this function searches for the header file and plots the data according to the header file of the data


        Para = np.loadtxt(params_fname[0],dtype='float',skiprows=1)
        dim = Para.shape #Dimensions of the loaded Array
        Winkeldata = np.array(Para[:,dim[1]-1])

        fig, axs = plt.subplots(2,3,gridspec_kw={'hspace': 0.3, 'wspace': 0.3})
        marker_size = 2
        symbol = '-o'

        #Slope and Offset are not repeating, so only once present in the file. Therefore no loop is needed 

        #Slope
        axs[1,2].plot(Winkeldata,np.array(Para[:,0]),'-o',markersize=marker_size)
        axs[1,2].set_title('Slope')

        #Offset
        axs[0,2].plot(Winkeldata,np.array(Para[:,1]),'-o',markersize=marker_size)
        axs[0,2].set_title('Offset')        

        try:
            if index_model == 2: #Lorentz
                num = int((dim[1]-3)/3) #Number of functions used for fitting
                for i_num in range(1,num+1):
                    #loop for the repeating parameters
                
                    #Linewidth
                    l4, = axs[0,1].plot(Winkeldata,np.array(Para[:,3*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[0,1].set_title('Linewidth [T]')

                    #Resonance Field
                    l5, = axs[1,0].plot(Winkeldata,np.array(Para[:,1+3*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[1,0].set_title('Resonance Field [T]')

                    #Amplitude
                    l6, = axs[1,1].plot(Winkeldata,np.array(Para[:,2+3*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[1,1].set_title('Amplitude [Arb. Units]')   

            elif index_model == 3: #Dyson
                num = int((dim[1]-3)/4) #Number of functions used for fitting
                for i_num in range(1,num+1):
                    #loop for the repeating parameters

                    #Alpha
                    axs[0,0].plot(Winkeldata,np.array(Para[:,4*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[0,0].set_title('Alpha')

                    #Linewidth
                    axs[0,1].plot(Winkeldata,np.array(Para[:,1+4*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[0,1].set_title('Linewidth [T]')

                    #Resonance Field
                    axs[1,0].plot(Winkeldata,np.array(Para[:,2+4*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[1,0].set_title('Resonance Field [T]')

                    #Amplitude
                    axs[1,1].plot(Winkeldata,np.array(Para[:,3+4*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[1,1].set_title('Amplitude [Arb. Units]')
        except Exception as e:
            print(e)
        plt.show()

    def plot_eigen(self,params_fname,index_model):
        Para = np.array(params_fname)
        dim = Para.shape #Dimensions of the loaded Array
        Winkeldata = np.array(Para[:,dim[1]-1])

        fig, axs = plt.subplots(2,3,gridspec_kw={'hspace': 0.3, 'wspace': 0.3})
        marker_size = 2
        symbol = '-o'

        #Slope and Offset are not repeating, so only once present in the file. Therefore no loop is needed 

        #Slope
        l1, = axs[1,2].plot(Winkeldata,np.array(Para[:,0]),'-o',markersize=marker_size)
        axs[1,2].set_title('Slope')

        #Offset
        l2, = axs[0,2].plot(Winkeldata,np.array(Para[:,1]),'-o',markersize=marker_size)
        axs[0,2].set_title('Offset')        
        try:
            if index_model == 2: #Lorentz
                num = int((dim[1]-3)/3) #Number of functions used for fitting
                for i_num in range(1,num+1):
                    #loop for the repeating parameters
                
                    #Linewidth
                    l4, = axs[0,1].plot(Winkeldata,np.array(Para[:,3*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[0,1].set_title('Linewidth [T]')

                    #Resonance Field
                    l5, = axs[1,0].plot(Winkeldata,np.array(Para[:,1+3*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[1,0].set_title('Resonanzfield [T]')

                    #Amplitude
                    l6, = axs[1,1].plot(Winkeldata,np.array(Para[:,2+3*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[1,1].set_title('Amplitude [Arb. Units]')                

            elif index_model == 3: #Dyson
                num = int((dim[1]-3)/4) #Number of functions used for fitting
                for i_num in range(1,num+1):
                    #loop for the repeating parameters

                    #Alpha
                    l3, = axs[0,0].plot(Winkeldata,np.array(Para[:,4*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[0,0].set_title('Alpha')

                    #Linewidth
                    l4, = axs[0,1].plot(Winkeldata,np.array(Para[:,1+4*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[0,1].set_title('Linewidth [T]')

                    #Resonance Field
                    l5, = axs[1,0].plot(Winkeldata,np.array(Para[:,2+4*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[1,0].set_title('Resonanzfield [T]')

                    #Amplitude
                    l6, = axs[1,1].plot(Winkeldata,np.array(Para[:,3+4*(i_num-1)+2]),'-o',markersize=marker_size)
                    axs[1,1].set_title('Amplitude [Arb. Units]')
        except Exception as e:
            print(e)
        plt.show()
