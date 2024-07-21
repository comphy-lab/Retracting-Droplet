import numpy as np
import os
import subprocess as sp
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import sys
import math
import pandas as pd

data = pd.read_csv('exp_data.csv')
t_exp = data.iloc[:,0]
theta_exp = data.iloc[:,1]

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True

def gettingFacets(filename):
    exe = ["./getFacet", filename]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    interface = []
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
                    midpoint = np.array([-(r1 + r2) / 2, (z1 + z2) / 2])
                 
                    interface.append(midpoint)
                    
                   
                    r0 = next((r for r, z in interface if z < 1e-1), None)

                    if r0 is not None:
                        # Sort the interface points based on their distance from (r0, 0)
                        interface.sort(key=lambda point: math.sqrt((point[0]-r0)**2 + point[1]**2))
                    # interface = sorted(interface, key=lambda x: x[0])
                    skip = True
    return interface
#---------------------------------------------------------------------------------------------------------------
def calculate_distances(interface):
    if not interface:
        return [], []
    
    interface_np = np.array(interface)
    diffs = np.diff(interface_np, axis=0)
    distances = np.sqrt(np.sum(diffs**2, axis=1))
    s = np.concatenate([[0], np.cumsum(distances)])
    
    r = interface_np[:, 0]
    z = interface_np[:, 1]
    
    # dr = np.diff(r)
    # dz = np.diff(z)
    dr_ds = np.gradient(r, s)
    dz_ds = np.gradient(z, s)
    
    # Calculate theta
    theta = np.arctan2(dz_ds, dr_ds)

    return s, theta
#---------------------------------------------------------------------------------------------------------------
nGFS=1000

for ti in range(nGFS):
    t = ti*(1e-1)
place = "test_young_50/snapshot-%5.4f" % 0.0000
place1 = "test_young_50/snapshot-%5.4f" % 15.0000
place2 = "test_young_50/snapshot-%5.4f" % 30.0000
place3 = "test_young_50/snapshot-%5.4f" % 45.0000
place4 = "test_young_50/snapshot-%5.4f" % 60.0000
interface = gettingFacets(place)
z = [point[1] for point in interface]
r = [point[0] for point in interface]
s, theta = calculate_distances(interface)

for i in range(len(theta)):
    # if theta[i] < 0:
    #     theta[i] += np.pi
    theta[i] = np.degrees(theta[i])   

interface1 = gettingFacets(place1)
z1 = [point[1] for point in interface1]
r1 = [point[0] for point in interface1]
s1, theta1 = calculate_distances(interface1)

for i in range(len(theta1)):
    # if theta1[i] < 0:
    #     theta1[i] += np.pi
    theta1[i] = np.degrees(theta1[i])       

interface2 = gettingFacets(place2)
z2 = [point[1] for point in interface2]
r2 = [point[0] for point in interface2]
s2, theta2 = calculate_distances(interface2)

for i in range(len(theta2)):
    # if theta2[i] < 0:
    #     theta2[i] += np.pi
    theta2[i] = np.degrees(theta2[i])   

interface3 = gettingFacets(place3)
z3 = [point[1] for point in interface3]
r3 = [point[0] for point in interface3]
s3, theta3 = calculate_distances(interface3)

for i in range(len(theta3)):
    # if theta3[i] < 0:
    #     theta3[i] += np.pi
    theta3[i] = np.degrees(theta3[i])   

interface4 = gettingFacets(place4)
z4 = [point[1] for point in interface4]
r4 = [point[0] for point in interface4]
s4, theta4 = calculate_distances(interface4)

for i in range(len(theta4)):
    # if theta4[i] < 0:
    #     theta4[i] += np.pi
    theta4[i] = np.degrees(theta4[i])           

fig, axs = plt.subplots(2, 2, figsize=(10, 6))
axs[0,0].plot(s, theta, linestyle='-')
axs[0,0].plot(s1, theta1, linestyle='-')
axs[0,0].plot(s2, theta2, linestyle='-')
axs[0,0].plot(s3, theta3, linestyle='-')
axs[0,0].plot(s4, theta4, linestyle='-')
axs[0,0].axvline(x = np.arcsin(0.3), color='r', linestyle='--')
# axs[0,0].set_ylim(0, 90)
# axs[0,0].set_xlim(0.,np.pi/2)
axs[0,0].set_xlabel('s', fontsize=20)
axs[0,0].set_ylabel('$\\theta$', fontsize=20)
# axs[0,0].set_xticks(fontsize=20)
# axs[0,0].set_yticks(fontsize=20)
axs[0,0].legend(['t=0', 't=15', 't=30', 't=60', 't=90'])
axs[0,0].grid(True)

axs[0,1].plot(z, s, linestyle='-')
axs[0,1].plot(z1, s1, linestyle='-')
axs[0,1].plot(z2, s2, linestyle='-')
axs[0,1].plot(z3, s3, linestyle='-')
axs[0,1].plot(z4, s4, linestyle='-') 
axs[0,1].axvline(x = (0.3), color='r', linestyle='--')
# axs[0,1].set_ylim(0, np.pi/2)
# axs[0,1].set_xlim(0.,0.9999)
axs[0,1].set_xlabel('h', fontsize=20)
axs[0,1].set_ylabel('s', fontsize=20)
# axs[0,1].set_xticks(fontsize=20)
# axs[0,1].set_yticks(fontsize=20)
axs[0,1].legend(['t=0', 't=15', 't=30', 't=60', 't=90'])
axs[0,1].grid(True)
    
axs[1,0].plot(z, theta, linestyle='-')
axs[1,0].plot(z1, theta1, linestyle='-')
axs[1,0].plot(z2, theta2, linestyle='-')
axs[1,0].plot(z3, theta3, linestyle='-')
axs[1,0].plot(z4, theta4, linestyle='-')   
axs[1,0].axvline(x = (0.3), color='r', linestyle='--')
# axs[1,0].set_ylim(0, 90)
# axs[1,0].set_xlim(0.,0.9999)
axs[1,0].set_xlabel('h', fontsize=20)
axs[1,0].set_ylabel('$\\theta$', fontsize=20)
# axs[1,0].set_xticks(fontsize=20)
# axs[1,0].set_yticks(fontsize=20)
axs[1,0].legend(['t=0', 't=15', 't=30', 't=60', 't=90'])
axs[1,0].grid(True)

axs[1,1].plot(r, z, linestyle='-')
axs[1,1].plot(r1, z1, linestyle='-')
axs[1,1].plot(r2, z2, linestyle='-')
axs[1,1].plot(r3, z3, linestyle='-')
axs[1,1].plot(r4, z4, linestyle='-')
axs[1,1].axhline(y = 0.3, color='r', linestyle='--')
axs[1,1].set_xlabel('r', fontsize=20)
axs[1,1].set_ylabel('h', fontsize=20)
axs[1,1].grid(True)

plt.savefig('no_young_test.png')
plt.show()



delZ = 2e-2
nGFS=1000
contact_angle,t_values = [], []

for ti in range(nGFS):
    t = ti*(1e-1)
    place = "test_young_50/snapshot-%5.4f" % t
    interface = gettingFacets(place)
    z = [point[1] for point in interface]
    r = [point[0] for point in interface]
    s, theta = calculate_distances(interface)

    first_theta = None
    for i in range(len(theta)):
        # if theta[i] < 0:
        #     theta[i] += np.pi
        theta[i] = np.degrees(theta[i])
        if abs(z[i] - 0.3) < delZ and first_theta is None:  
            contact_angle.append(theta[i])
            t_values.append(t)
            found = True
            break

# contact_angle1, t_values1 = [], []

# for ti in range(nGFS):
#     t = ti*(1e-1)
#     place = "test_25/snapshot-%5.4f" % t
#     interface = gettingFacets(place)
#     z = [point[1] for point in interface]
#     r = [point[0] for point in interface]
#     s, theta = calculate_distances(interface)

#     first_theta = None
#     for i in range(len(theta)):
#         if theta[i] < 0:
#             theta[i] += np.pi
#         theta[i] = np.degrees(theta[i])
#         if abs(z[i] - 0.3) < delZ and first_theta is None:  
#             contact_angle1.append(theta[i])
#             t_values1.append(t)
#             found = True
#             break


t0 = 3.38
t_values = [i * t0 for i in t_values]
# t_values1 = [i * t0 for i in t_values1]
# print(contact_angle)

plt.figure()
plt.plot(t_values, contact_angle, marker='o', linestyle='-')
# plt.plot(t_values1, contact_angle1, marker='s', linestyle='-')
# plt.plot(t_exp, theta_exp, marker = '*')
# plt.ylim(0, 100)
# plt.xlim(2,18)
plt.grid(True)
plt.xlabel('t (ms)', fontsize=18)
plt.ylabel('$\\theta$', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
# plt.title('$\\theta_{s} = 30 ^\circ$', fontsize=20)
plt.legend(['Hysteresis model', 'Experimental Data'])
plt.subplots_adjust(bottom=0.15)
plt.savefig('contact_angle_young_70.png')
plt.show()