#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Processing for UNH-RVAT 3D OpenFOAM simulation.

by Pete Bachant (petebachant@gmail.com)

"""
from __future__ import division, print_function
import matplotlib.pyplot as plt
import re
import numpy as np
import os
from styleplot import setpltparams
import sys
import foampy
import fdiff
    
exp_path = "/media/pete/External 2/Research/Experiments/2014 Spring RVAT Re dep"

#setpltparams()

# Some constants
R = 0.13
U = 1.0
H = 0.26
D = 2*R
A = np.pi*R**2
rho = 1000.0

ylabels = {"meanu" : r"$U/U_\infty$",
           "stdu" : r"$\sigma_u/U_\infty$",
           "meanv" : r"$V/U_\infty$",
           "meanw" : r"$W/U_\infty$",
           "meanuv" : r"$\overline{u'v'}/U_\infty^2$"}
    
def loadwake(time):
    """Loads wake data and returns y/R and statistics."""
    # Figure out if time is an int or float
    if not isinstance(time, str):
        if time % 1 == 0:
            folder = str(int(time))
        else:
            folder = str(time)
    else:
        folder = time
    flist = os.listdir("postProcessing/sets/"+folder)
    data = {}
    for fname in flist:
        fpath = "postProcessing/sets/"+folder+"/"+fname
        z_H = float(fname.split("_")[1])/R
        data[z_H] = np.loadtxt(fpath, unpack=True)
    return data
    
def calcwake(t1=0.0):
    times = os.listdir("postProcessing/sets")
    times = [float(time) for time in times]
    times.sort()
    times = np.asarray(times)
    data = loadwake(times[0])
    z_H = np.asarray(sorted(data.keys()))
    y_R = data[z_H[0]][0]/R
    # Find first timestep from which to average over
    i = np.where(times==t1)[0][0]
    t = times[i:]
    # Assemble 3-D arrays, with time as first index
    u = np.zeros((len(t), len(z_H), len(y_R)))
    v = np.zeros((len(t), len(z_H), len(y_R)))
    w = np.zeros((len(t), len(z_H), len(y_R)))
    xvorticity = np.zeros((len(t), len(z_H), len(y_R)))
    # Loop through all times
    for n in range(len(t)):
        data = loadwake(t[n])
        for m in range(len(z_H)):
            u[n,m,:] = data[z_H[m]][1]
            v[n,m,:] = data[z_H[m]][2]
            w[n,m,:] = data[z_H[m]][3]
            xvorticity[n,m,:] = data[z_H[m]][4]
    meanu = u.mean(axis=0)
    meanv = v.mean(axis=0)
    meanw = w.mean(axis=0)
    xvorticity = xvorticity.mean(axis=0)
    return {"meanu" : meanu,
            "meanv" : meanv,
            "meanw" : meanw,
            "xvorticity" : xvorticity,
            "y/R" : y_R, 
            "z/H" : z_H}
        
    
def plotwake(plotlist=["meanu"], t1=3.0, save=False, savepath="", 
             savetype=".pdf"):
    data = calcwake(t1=t1)
    y_R = data["y/R"]
    z_H = data["z/H"]
    u = data["meanu"]
    v = data["meanv"]
    w = data["meanw"]
    xvorticity = data["xvorticity"]
    def turb_lines(half=False):
        if half:
            plt.hlines(0.5, -1, 1, linestyles='solid', linewidth=2)
            plt.vlines(-1, 0, 0.5, linestyles='solid', linewidth=2)
            plt.vlines(1, 0, 0.5, linestyles='solid', linewidth=2)
        else:
            plt.hlines(0.5, -1, 1, linestyles='solid', colors='gray',
                       linewidth=3)
            plt.hlines(-0.5, -1, 1, linestyles='solid', colors='gray',
                       linewidth=3)
            plt.vlines(-1, -0.5, 0.5, linestyles='solid', colors='gray',
                       linewidth=3)
            plt.vlines(1, -0.5, 0.5, linestyles='solid', colors='gray',
                       linewidth=3)
    if "meanu" in plotlist or "all" in plotlist:
        plt.figure(figsize=(9,8))
        cs = plt.contourf(y_R, z_H, u, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.12)
        cb.set_label(r'$U/U_{\infty}$')
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(1)
        plt.tight_layout()
    if "meanv" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y/0.5, z, v, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        plt.tight_layout()
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.22)
        cb.set_label(r'$V/U_{\infty}$')
        #turb_lines()
        ax = plt.axes()
        ax.set_aspect(1)
        plt.grid(True)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
    if "v-wquiver" in plotlist or "all" in plotlist:
        # Make quiver plot of v and w velocities
        plt.figure(figsize=(10,5))
        Q = plt.quiver(y_R, z_H, v, w, angles='xy')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        plt.quiverkey(Q, 0.75, 0.2, 0.1, r'$0.1$ m/s',
                   labelpos='E',
                   coordinates='figure',
                   fontproperties={'size': 'small'})
        plt.tight_layout()
        plt.hlines(0.5, -1, 1, linestyles='solid', colors='r',
                   linewidth=2)
        plt.vlines(-1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=2)
        plt.vlines(1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=2)
        ax = plt.axes()
        ax.set_aspect(1)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'v-wquiver'+savetype)
    if "xvorticity" in plotlist or "all" in plotlist:
        plt.figure(figsize=(6, 8))
        cs = plt.contourf(y_R, z_H, xvorticity, 10, cmap=plt.cm.coolwarm)
                          #levels=np.linspace(-2.5,2.5,21))
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/R$')
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.08)
#        cb.set_ticks(np.linspace(-2.5,2.5,11), update_ticks=True)
        cb.set_label(r"$\Omega_x$")
        #turb_lines()
        ax = plt.axes()
        ax.set_aspect(1)
        plt.tight_layout()
        if save:
            plt.savefig(savepath+'/xvorticity_miniHAWT'+savetype)
    if "meancomboquiv" in plotlist or "all" in plotlist:
        plt.figure(figsize=(6, 8))
        # Add contours of mean velocity
        cs = plt.contourf(y_R, z_H, u, 20, cmap=plt.cm.coolwarm)
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.08)
        cb.set_label(r'$U/U_{\infty}$')
        plt.hold(True)
        # Make quiver plot of v and w velocities
        Q = plt.quiver(y_R, z_H, v, w, angles='xy', width=0.0022)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/R$')
        plt.quiverkey(Q, 0.8, 0.19, 0.1, r'$0.1 U_\infty$',
                      labelpos='E',
                      coordinates='figure',
                      fontproperties={'size': 'small'})
        ax = plt.axes()
        ax.set_aspect(1)
        plt.tight_layout()
        if save:
            plt.savefig(savepath+"\\meancomboquiv_miniHAWT"+savetype)
    plt.show()
        
def plotexpwake(Re_D, quantity, z_H=0.0, save=False, savepath="", 
                savetype=".pdf", newfig=True, marker="--ok",
                fill="none", figsize=(10, 5)):
    """Plots the transverse wake profile of some quantity. These can be
      * meanu
      * meanv
      * meanw
      * stdu
    """
    U = Re_D/1e6
    label = "Exp."
    folder = exp_path + "/Wake/U_" + str(U) + "/Processed/"
    z_H_arr = np.load(folder + "z_H.npy")
    i = np.where(z_H_arr==z_H)
    q = np.load(folder + quantity + ".npy")[i]
    y_R = np.load(folder + "y_R.npy")[i]
    if newfig:
        plt.figure(figsize=figsize)
    plt.plot(y_R, q/U, marker, markerfacecolor=fill, label=label)
    plt.xlabel(r"$y/R$")
    plt.ylabel(ylabels[quantity])
    plt.grid(True)
    plt.tight_layout()

def main():
    p = "Google Drive/temp"
    if "linux" in sys.platform:
        p = "/home/pete/" + p
    elif "win" in sys.platform:
        p = "C:/Users/Pete/" + p
    plt.close("all")
    
    plotwake(plotlist=["xvorticity", "meancomboquiv"], t1=2.0, 
             save=True, savepath=p)
#    calcwake()

if __name__ == "__main__":
    main()
