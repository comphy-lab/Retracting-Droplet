# Author: Vatsal Sanjay
# Edited by: Aman Bhargava
# vatsalsanjay@gmail.com, amanbhargava2000@gmail.com
# Physics of Fluids
# Last updated: 02-Feb-2024

import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import sys

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True

def gettingFacets(filename):
    exe = ["./getFacet", filename]
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
                    segs.append(((r1, z1),(r2,z2)))
                    segs.append(((-r1, z1),(-r2,z2)))
                    skip = True
    return segs

def gettingfield(filename):
    exe = ["./getData", filename, str(zmin), str(rmin), str(zmax), str(rmax), str(nr)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    # print(temp2) #debugging
    Rtemp, Ztemp, D2temp, veltemp, momtemp = [],[],[],[],[]
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
    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    D2 = np.asarray(D2temp)
    vel = np.asarray(veltemp)
    mom = np.asarray(momtemp)
    nz = int(len(Z)/nr)
    # print("nr is %d %d" % (nr, len(R))) # debugging
    print("nz is %d" % nz)
    R.resize((nz, nr))
    Z.resize((nz, nr))
    D2.resize((nz, nr))
    vel.resize((nz, nr))
    mom.resize((nz,nr))

    return R, Z, D2, vel, mom, nz
# ----------------------------------------------------------------------------------------------------------------------


Ldomain = 9.0
nr = 256
rmin, rmax, zmin, zmax = [0.0, Ldomain, 0.0, Ldomain]
lw = 4

folder = 'hysteresis_test_fixed'  # output folder

if not os.path.isdir(folder):
    os.makedirs(folder)

nGFS=10000

for ti in range(nGFS):
    t = ti*(1e-2)
    print("Time is %f" % t)

    place = "test_hysteresis_fixed/snapshot-%5.4f" % t
    name = "%s/%8.8d.png" %(folder, int(t*1000000))

    if not os.path.exists(place):
        print("%s File not found!" % place)
    else:
        if os.path.exists(name):
            print("%s Image present!" % name)
        else:

            segs = gettingFacets(place)
            if (len(segs) == 0):
                print("Problem in the available file %s" % place)
            else:
                R, Z, D2, vel, mom, nz = gettingfield(place)
                zminp, zmaxp, rminp, rmaxp = Z.min(), Z.max(), R.min(), R.max()
                # print(zminp, zmaxp, rminp, rmaxp)
                # Part to plot
                AxesLabel, TickLabel = [50, 20]
                fig, ax = plt.subplots()
                fig.set_size_inches(19.20, 10.80)

                ## Drawing Facets
                line_segments = LineCollection(segs, linewidths=4, colors='#1a9850', linestyle='solid')
                ax.add_collection(line_segments)

                ## D
                cntrl1 = ax.imshow(D2, cmap="hot_r", interpolation='None', origin='lower', extent=[rminp, rmaxp, zminp, zmaxp], vmax = -1.0, vmin = -3.0)
                ## V
                cntrl2 = ax.imshow(mom, interpolation='None', cmap="Blues", origin='lower', extent=[-rminp, -rmaxp, zminp, zmaxp], vmax = 1e0, vmin = 0.)

                ax.plot([0, 0], [zmin, zmax],'-.',color='grey',linewidth=lw)
                # ax.plot([rmin, rmax], [0, 0],'-',color='grey',linewidth=lw/2)
                ax.plot([-rmax, -rmax], [zmin, zmax],'-',color='black',linewidth=lw)
                ax.plot([-rmax, rmax], [zmin, zmin],'-',color='black',linewidth=lw)
                ax.plot([-rmax, rmax], [zmax, zmax],'-',color='black',linewidth=lw)
                ax.plot([rmax, rmax], [zmin, zmax],'-',color='black',linewidth=lw)
                # add line where experiments are being measured
                plt.axhline(y=0.3, color='k', linestyle='--')

                ax.set_aspect('equal')
                ax.set_xlim(-rmax, rmax)
                ax.set_ylim(zmin, zmax)
                ax.set_title('$t/t_\gamma$ = %5.4f' % t, fontsize=TickLabel)

                l, b, w, h = ax.get_position().bounds
                cb1 = fig.add_axes([l+0.55*w, b-0.05, 0.40*w, 0.03])
                c1 = plt.colorbar(cntrl1,cax=cb1,orientation='horizontal')
                c1.set_label('$\log_{10}\left(\\varepsilon_\eta\\right)$',fontsize=TickLabel, labelpad=5)
                c1.ax.tick_params(labelsize=TickLabel)
                c1.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
                cb2 = fig.add_axes([l+0.05*w, b-0.05, 0.40*w, 0.03])
                c2 = plt.colorbar(cntrl2,cax=cb2,orientation='horizontal')
                c2.ax.tick_params(labelsize=TickLabel)
                #c2.set_label('$\|v_i\|/V_\gamma$',fontsize=TickLabel)
                c2.set_label('$\\rho U$',fontsize=TickLabel)
                c2.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
                ax.axis('off')
                # plt.show()
                plt.savefig(name, bbox_inches="tight")
                plt.close()
