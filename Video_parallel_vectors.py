# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last updated: Oct 27, 2024

import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import multiprocessing as mp
from functools import partial
import argparse  # Add at top with other imports

import matplotlib.colors as mcolors
# custom_colors = ["white", "#DA8A67", "#A0522D", "#400000"]
# custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_hot", custom_colors)

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

def gettingFacets(filename,includeCoat='true'):
    exe = ["./getFacet", filename, includeCoat]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    segs.append(((r1, z1),(r2, z2)))
                    segs.append(((-r1, z1),(-r2, z2)))
                    skip = True
    return segs

def gettingfield(filename, zmin, rmin, zmax, rmax, nr, t):
    exe = ["./getData_0.0075", filename, str(zmin), str(rmin), str(zmax), str(rmax), str(nr)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    # print(temp2) #debugging
    Rtemp, Ztemp, D2temp, veltemp, momtemp, momxtemp, momytemp  = [],[],[],[],[],[],[]

    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Ztemp.append(float(temp3[0]))
            Rtemp.append(float(temp3[1]))
            D2temp.append(float(temp3[2]))
            veltemp.append(float(temp3[3]))
            momtemp.append(float(temp3[4]))
            momxtemp.append(float(temp3[5]))
            momytemp.append(float(temp3[6]))

    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    D2 = np.asarray(D2temp)
    vel = np.asarray(veltemp)
    mom = np.asarray(momtemp)
    momx = np.asarray(momxtemp)
    momy = np.asarray(momytemp)
    nz = int(len(Z)/nr)

    # print("nr is %d %d" % (nr, len(R))) # debugging
    print("nz is %d" % nz)
    print("Time is %f" % t)

    R.resize((nz, nr))
    Z.resize((nz, nr))
    D2.resize((nz, nr))
    vel.resize((nz, nr))
    mom.resize((nz,nr))
    momx.resize((nz,nr))
    momy.resize((nz,nr))

    return R, Z, D2, vel, mom, momx, momy, nz
# ----------------------------------------------------------------------------------------------------------------------

def process_timestep(ti, folder, nGFS, GridsPerR, rmin, rmax, zmin, zmax, lw):
    t = ti*(1e-1)
    nr = 256
    place = "999/snapshot-%5.4f" % t
    name = "%s/%8.8d.pdf" %(folder, int(t*1000000))
    if not os.path.exists(place):
        print(f"{place} File not found!")
        return

    if os.path.exists(name):
        print(f"{name} Image present!")
        return

    segs = gettingFacets(place)
    if (len(segs) == 0):
        print("Problem in the available file %s" % place)
    else:
        R, Z, D2, vel, mom, momx, momy, nz = gettingfield(place, zmin, rmin, zmax, rmax, nr, t)
        zminp, zmaxp, rminp, rmaxp = Z.min(), Z.max(), R.min(), R.max()

        # Plotting
        AxesLabel, TickLabel = 50, 20
        fig, ax = plt.subplots()
        fig.set_size_inches(19.20, 10.80)

        rmax = 7.0
        zmax = 3.0

        ax.plot([0, 0], [zmin, zmax], '-.', color='grey', linewidth=lw)
        ax.plot([-rmax, -rmax], [zmin, zmax], '-', color='black', linewidth=lw)
        ax.plot([-rmax, rmax], [zmin, zmin], '-', color='black', linewidth=lw)
        ax.plot([-rmax, rmax], [zmax, zmax], '-', color='black', linewidth=lw)
        ax.plot([rmax, rmax], [zmin, zmax], '-', color='black', linewidth=lw)
        # add line where experiments are being measured
        plt.axhline(y=0.3, color='red', linestyle='--')

        line_segments = LineCollection(segs, linewidths=4, colors='black', linestyle='solid')
        ax.add_collection(line_segments)

        ## Dissipation
        cntrl1 = ax.imshow(D2, cmap="hot_r", interpolation='None', origin='lower', extent=[rminp, rmaxp, zminp, zmaxp], vmax = -1.0, vmin = -3.0)
        ## Momentum
        cntrl2 = ax.imshow(mom, interpolation='None', cmap="Blues", origin='lower', extent=[-rminp, -rmaxp, zminp, zmaxp], vmax = 1e0, vmin = 0.)
        
        n = 3  # Adjust this value to control the density
        R_reduced = R[::n, ::n]
        Z_reduced = Z[::n, ::n]
        momx_reduced = momx[::n, ::n]
        momy_reduced = momy[::n, ::n]

        magnitude = np.sqrt(momx_reduced**2 + momy_reduced**2)

        mask = magnitude > 0.01

        R_filtered = R_reduced[mask]
        Z_filtered = Z_reduced[mask]
        momx_filtered = momx_reduced[mask]
        momy_filtered = momy_reduced[mask]

        # Plot vectors
        ax.quiver(-R_filtered, Z_filtered, -momy_filtered, momx_filtered, color='#1F77B4', scale=10, headwidth=3, headlength=4, headaxislength=4, width=0.002)
        ax.quiver(R_filtered, Z_filtered, momy_filtered, momx_filtered, color='#1F77B4', scale=10, headwidth=3, headlength=4, headaxislength=4, width=0.002)     
        # ax.quiver(R, Z, ux, uy, color='black', scale=1000)

        ax.set_aspect('equal')
        ax.set_xlim(-rmax, rmax)
        ax.set_ylim(0, zmax)
        ax.set_title('$t/t_\\gamma$ = %5.4f' % t, fontsize=TickLabel)

        l, b, w, h = ax.get_position().bounds
        # cb1 = fig.add_axes([l+0.55*w, b-0.05, 0.40*w, 0.03])
        # c1 = plt.colorbar(cntrl1,cax=cb1,orientation='horizontal')
        # c1.set_label('$\\log_{10}\\left(\\varepsilon_\\eta\\right)$',fontsize=36, labelpad=5)
        # c1.ax.tick_params(labelsize=TickLabel)
        # c1.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
        # cb2 = fig.add_axes([l+0.05*w, b-0.05, 0.40*w, 0.03])
        # c2 = plt.colorbar(cntrl2,cax=cb2,orientation='horizontal')
        # c2.ax.tick_params(labelsize=TickLabel)
        # #c2.set_label('$\|v_i\|/V_\gamma$',fontsize=TickLabel)
        # c2.set_label('$\\rho U$',fontsize=36)
        # c2.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
        ax.axis('off')
        # plt.show()
        plt.savefig(name, bbox_inches="tight")
        plt.close()

def main():
    LDomain = 16.0
    # Get number of CPUs from command line argument, or use all available
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--CPUs', type=int, default=mp.cpu_count(), help='Number of CPUs to use')
    parser.add_argument('--nGFS', type=int, default=2500, help='Number of restart files to process')
    parser.add_argument('--ZMAX', type=float, default=LDomain, help='Maximum Z value')
    parser.add_argument('--RMAX', type=float, default=LDomain, help='Maximum R value')
    parser.add_argument('--ZMIN', type=float, default=0.0, help='Minimum Z value')
    parser.add_argument('--RMIN', type=float, default=0.0, help='Minimum R value')
    args = parser.parse_args()

    CPUStoUse = args.CPUs
    nGFS = args.nGFS
    ZMAX = args.ZMAX
    RMAX = args.RMAX
    ZMIN = args.ZMIN
    RMIN = args.RMIN
    
    num_processes = CPUStoUse
    rmin, rmax, zmin, zmax = [RMIN, RMAX, ZMIN, ZMAX]
    GridsPerR = 128


    lw = 4
    folder = '999-vectors1'

    if not os.path.isdir(folder):
        os.makedirs(folder)

    # Create a pool of worker processes
    with mp.Pool(processes=num_processes) as pool:
        # Create partial function with fixed arguments
        process_func = partial(process_timestep, 
                             folder=folder, nGFS=nGFS,
                             GridsPerR=GridsPerR, rmin=rmin, rmax=rmax, 
                             zmin=zmin, zmax=zmax, lw=lw)
        # Map the process_func to all timesteps
        pool.map(process_func, range(nGFS))

if __name__ == "__main__":
    main()
