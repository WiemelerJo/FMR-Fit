import matplotlib as mp
mp.use('QT5Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
import math as m

class Params_Plot(object):
    """docstring for Params_Plot"""
    def __init__(self, params_fname,BOOL):
        #QThread.__init__(self)
        super().__init__()
        #self.params_fname = params_fname
        if BOOL == 'eigen':
            self.plot_eigen(params_fname)
        else:
            self.plot(params_fname)

    def plot(self,params_fname):
        Para = np.loadtxt(params_fname[0],dtype='float',skiprows=1)
        dim = Para.shape #Dimensions of the loaded Array
        num = int((dim[1]-1)/6) #Number of functions used for fitting
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

        for i_num in range(0,num+1):
            print(i_num+1)
            #loop for the repeating parameters

            #Alpha
            axs[0,0].plot(Winkeldata,np.array(Para[:,4*i_num-1+2]),'-o',markersize=marker_size)
            axs[0,0].set_title('Alpha'+str(i_num+1))

            #Linewidth
            axs[0,1].plot(Winkeldata,np.array(Para[:,1+4*i_num-1+2]),'-o',markersize=marker_size)
            axs[0,1].set_title('Linewidth [T]'+str(i_num+1))

            #Resonance Field
            axs[1,0].plot(Winkeldata,np.array(Para[:,2+4*i_num-1+2]),'-o',markersize=marker_size)
            axs[1,0].set_title('Resonanzfield [T]'+str(i_num+1))

            #Amplitude
            axs[1,1].plot(Winkeldata,np.array(Para[:,3+4*i_num-1+2]),'-o',markersize=marker_size)
            axs[1,1].set_title('Amplitude [Arb. Units]'+str(i_num+1))
        plt.show()

    def plot_eigen(self,params_fname):
        Para = np.array(params_fname)
        dim = Para.shape #Dimensions of the loaded Array
        num = int((dim[1]-1)/6) #Number of functions used for fitting
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

        for i_num in range(0,num+1):
            print(i_num+1)
            #loop for the repeating parameters

            #Alpha
            axs[0,0].plot(Winkeldata,np.array(Para[:,4*i_num-1+2]),'-o',markersize=marker_size)
            axs[0,0].set_title('Alpha'+str(i_num+1))

            #Linewidth
            axs[0,1].plot(Winkeldata,np.array(Para[:,1+4*i_num-1+2]),'-o',markersize=marker_size)
            axs[0,1].set_title('Linewidth [T]'+str(i_num+1))

            #Resonance Field
            axs[1,0].plot(Winkeldata,np.array(Para[:,2+4*i_num-1+2]),'-o',markersize=marker_size)
            axs[1,0].set_title('Resonanzfield [T]'+str(i_num+1))

            #Amplitude
            axs[1,1].plot(Winkeldata,np.array(Para[:,3+4*i_num-1+2]),'-o',markersize=marker_size)
            axs[1,1].set_title('Amplitude [Arb. Units]'+str(i_num+1))
        plt.show()
