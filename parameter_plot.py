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
        Alpha_param = np.array(Para[:,0])
        dB_param = np.array(Para[:,1])
        R_param = np.array(Para[:,2])
        A_param = np.array(Para[:,3])
        slope_param = np.array(Para[:,4])
        offset_param = np.array(Para[:,5])
        Winkeldata = np.array(Para[:,6])

        fig, axs = plt.subplots(2,3,gridspec_kw={'hspace': 0.3, 'wspace': 0.3})
        marker_size = 2
        symbol = '-o'

        axs[0,0].plot(Winkeldata,Alpha_param,'-o',markersize=marker_size)
        axs[0,0].set_title('Alpha')
        axs[0,1].plot(Winkeldata,dB_param,'-o',markersize=marker_size)
        axs[0,1].set_title('Linewidth [T]')
        axs[1,0].plot(Winkeldata,R_param,'-o',markersize=marker_size)
        axs[1,0].set_title('Resonanzfield [T]')
        axs[1,1].plot(Winkeldata,A_param,'-o',markersize=marker_size)
        axs[1,1].set_title('Amplitude [Arb. Units]')
        axs[1,2].plot(Winkeldata,slope_param,'-o',markersize=marker_size)
        axs[1,2].set_title('Slope')
        axs[0,2].plot(Winkeldata,offset_param,'-o',markersize=marker_size)
        axs[0,2].set_title('Offset')
        plt.show()

    def plot_eigen(self,params_fname):
        Para = np.array(params_fname)
        Alpha_param = np.array(Para[:,0])
        dB_param = np.array(Para[:,1])
        R_param = np.array(Para[:,2])
        A_param = np.array(Para[:,3])
        slope_param = np.array(Para[:,4])
        offset_param = np.array(Para[:,5])
        Winkeldata = np.array(Para[:,6])

        fig, axs = plt.subplots(2,3,gridspec_kw={'hspace': 0.3, 'wspace': 0.3})
        marker_size  =2
        symbol = '-o'

        axs[0,0].plot(Winkeldata,Alpha_param,'-o',markersize=marker_size)
        axs[0,0].set_title('Alpha')
        axs[0,1].plot(Winkeldata,dB_param,'-o',markersize=marker_size)
        axs[0,1].set_title('Linewidth [T]')
        axs[1,0].plot(Winkeldata,R_param,'-o',markersize=marker_size)
        axs[1,0].set_title('Resonanzfield [T]')
        axs[1,1].plot(Winkeldata,A_param,'-o',markersize=marker_size)
        axs[1,1].set_title('Amplitude [Arb. Units]')
        axs[1,2].plot(Winkeldata,slope_param,'-o',markersize=marker_size)
        axs[1,2].set_title('Slope')
        axs[0,2].plot(Winkeldata,offset_param,'-o',markersize=marker_size)
        axs[0,2].set_title('Offset')
        plt.show()
