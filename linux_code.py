from __future__ import division
import csv
from collections import defaultdict
import subprocess as sp
import os
import shutil
import sys
import string
import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

run = "xfoil"

def Xfoil(name, Re, Alpha , Mach ):
    def Cmd(cmd):
        ps.stdin.write(cmd+'\n')
    try:
        os.remove(name+'.log')
    except :
        pass
    
    #Run XFOIL
    
    ps = sp.Popen(run ,stdin=sp.PIPE,stderr=sp.PIPE,stdout=sp.PIPE, shell = True)
    ps.stderr.close()
    
    # XFOIL Commands
    
    Cmd('load '+name+'.csv')
    Cmd('OPER')
    Cmd('visc '+str(Re))
    Cmd('mach'+str(Mach))
    Cmd('alfa' + str(Alpha))
    Cmd('cpwr')
    Cmd(name+'CP'+'.csv')
    Cmd(' ')     
    Cmd('quit')  
    
    ps.stdout.close()
    ps.stdin.close()
    ps.wait()
        
    #Getting X and CP in Array from File Generated in XFOIL

    filename = name+'CP'+'.csv'
    f = open(filename, 'r')

    flines = f.readlines()
    xvalue=[]
    cpvalue=[]
    for i in range(1,len(flines)):
    	flines[i] = flines[i].strip()
        xvalue.append(flines[i].split(' ')[0])
        if(flines[i].split(" ")[2]):
        	cpvalue.append(flines[i].split(' ')[2])
        elif(flines[i].split(" ")[3]):
        	cpvalue.append(flines[i].split(' ')[3])
        else:
        	cpvalue.append(flines[i].split(' ')[4])
    cpvalue = map(float, cpvalue)
    xvalue = map(float, xvalue)
    cp_x = []
    cp_x.append(cpvalue)
    cp_x.append(xvalue)

    return cp_x

def ReadXY(name):

    #Getting X and Y in Array
    
    filenameXY = name+'.csv'
    xy = open(filenameXY, 'r')

    flinesxy = xy.readlines()
    xvalues = []
    yvalues = []
    for j in range(1,len(flinesxy)):
    	flinesxy[j] = flinesxy[j].strip()
    	xvalues.append(flinesxy[j].split(' ')[0])
        yvalues.append(flinesxy[j].split(' ')[1])
    xvalues = map(float, xvalues)
    yvalues = map(float, yvalues)
    xy_array = []
    xy_array.append(xvalues)
    xy_array.append(yvalues)
    return xy_array

def create_xy(x_y_array):
    
    #Saving Modified File
    
 	ascii.write([x_y_array[0],x_y_array[1]], 'Modified.csv', names=['Modified', 'y'])



########################################################################################################################################

def algorithm(target_cp, current_cp, airfoil_xy, mach):
#airfoil_xy is a 2d array- x and y
#It is assumed that the length of the target/current cp are same and the no of x coordinates is also the same
#ambient conditions
    p_0 = 1.0125 * 10**5
    t_0 = 288
    gam = 1.4
    rho = 1.225
    speed_sound = (gam * t_0 * 287)**0.5
    q_0 = 0.5 * rho * (mach*speed_sound)**2 
    ep = 0.15 # for subsonic flow
    dt = 0.0001
#step 1: calculating pressures from cp
    current_press = []
    target_press = []
    for i in range(len(target_cp)):
        target_press.append(target_cp[i]*q_0 + p_0)
        current_press.append(current_cp[i]*q_0 + p_0)

    

#Mostly this is for the subsonic flow
#step2: calculating slopes of the normals
    slope = []
    for i in range(len(target_cp)-1):
        if airfoil_xy[0][i+1]-airfoil_xy[0][i] == 0:
            slope.append(1)
        else:
            slope.append((airfoil_xy[1][i+1]-airfoil_xy[1][i])/(airfoil_xy[0][i+1]-airfoil_xy[0][i]))
        
#think about the last point
#step3: calculating the virtual velocities
    virtual_y = []
    for i in range(len(target_cp)-1):
        diff = (target_press[i] - current_press[i])/rho
        if slope[i] == 0:
        	slope[i] = 1
        if i == 0:
            virtual_y.append(0)
        if i!=0 and diff <= 0: 
            virtual_y.append(abs((1/(slope[i]**2 + 1)) * diff)**0.5)        
        if i!=0 and diff>0:
            virtual_y.append(-abs(((1/(slope[i]**2+1)*diff))**0.5))


#calculating the relaxation factor
    omega = []
    for i in range(len(target_cp)-1):
        omega.append(ep*(1/speed_sound)*(abs((target_press[i] - current_press[i])/rho
)**0.5))

#calculating displacement parallel to the y axis
    dy = []
    for i in range(len(target_cp)-1):
        dy.append(omega[i]*virtual_y[i]*dt)
#Eliminating the smoothing part and i am calculating over leading and the trailing edges also-could be a problem

    for i in range(len(target_cp)-1):
        
        airfoil_xy[1][i] = airfoil_xy[1][i] + dy[i]
    for i in range(2, len(target_cp)-2):
        airfoil_xy[1][i] = (airfoil_xy[1][i-2]+airfoil_xy[1][i-1] + airfoil_xy[1][i] + airfoil_xy[1][i+1]+airfoil_xy[1][i+2])/5

    return(airfoil_xy)
#############################################################################################################################

def get_cp_on_x(cp_x, airfoil_x):

    #cp_x is a 2d array- consisting of cp and its x coordinates, whereas airfoil is a 1d array consisting of starting airfoil x coordinates
    
    tell = 0
    while cp_x[1][tell]!=0.00000:
        tell = tell+1

    upper_tar = [[0 for i in range(tell+1)] for y in range(2)]
    lower_tar = [[0 for i in range(len(cp_x[0])-tell)] for y in range(2)]
    

    for i in range(tell+1): 
        upper_tar[0][i] = cp_x[0][i]
        upper_tar[1][i] = cp_x[1][i]

    
    for i in range(len(cp_x[0])-tell):
        lower_tar[0][i] = cp_x[0][tell+i]
        lower_tar[1][i] = cp_x[1][tell+i]

    

    tel3 = 0
    while airfoil_x[tel3]!=0:
        tel3 = tel3+1
    airfoil_x_upper = [0 for i in range(tel3+1)]
    airfoil_x_lower = [0 for i in range(len(airfoil_x)-tel3)]
    for i in range(tel3+1): 
        airfoil_x_upper[i] = airfoil_x[i]
        

    
    for i in range(len(airfoil_x)-tel3):
        airfoil_x_lower[i] = airfoil_x[tel3+i]
        


    len_airfoil_lower = len(airfoil_x_lower)
    len_airfoil_upper = len(airfoil_x_upper)
    len_lower_x = len(lower_tar[1])
    len_upper_x = len(upper_tar[1])
    cp_airfoil_upper = []
    cp_airfoil_lower = []
    for i in range(len_airfoil_lower):
        for j in range(1,len_lower_x):
            if airfoil_x_lower[i] <= lower_tar[1][j] and airfoil_x_lower[i] >= lower_tar[1][j-1]:
                data = lower_tar[0][j] + (airfoil_x_lower[i]-lower_tar[1][j]) * (lower_tar[0][j]-lower_tar[0][j-1])/(lower_tar[1][j]-lower_tar[1][j-1])
                cp_airfoil_lower.append(data)
                break
    for i in range(len_airfoil_upper):
        for j in range(1,len_upper_x):
            if airfoil_x_upper[i] >= upper_tar[1][j] and airfoil_x_upper[i] <= upper_tar[1][j-1]:
                data = upper_tar[0][j] + (airfoil_x_upper[i]-upper_tar[1][j]) * (upper_tar[0][j]-upper_tar[0][j-1])/(upper_tar[1][j]-upper_tar[1][j-1])
                cp_airfoil_upper.append(data)
                break
    del cp_airfoil_upper[-1]

    final = cp_airfoil_upper + cp_airfoil_lower
    return final

target_cp_x = Xfoil('n0009', '5000000' , '0', '0.2')                       #Generating Target Cp
target_airfoil = ReadXY('n0009')                                           #Reading X-Y of Target Airfoil
initial_airfoil = ReadXY('n0015')                                          #Reading X-Y of Initial Airfoil
modified_x_y = ReadXY('Modified')                                          #Reading X-Y of Modified Airfoil
transfered_cp_x = get_cp_on_x(target_cp_x,modified_x_y[0])                 #Transfering X coordinates of Modified over Target
l2array = []
for x in xrange(1,170):   
    #Running Iteration for 170 times
    current_cp = Xfoil('Modified', '5000000', '0', '0.2')                          #Generating Modified Cp
    modified_x_y = algorithm(transfered_cp_x, current_cp[0], modified_x_y, 0.2)    #Running Alogithm on Modified Cp
    create_xy(modified_x_y)                                                        #Creating New Modified File
    target_l2_norm = Xfoil('n0009', '5000000' , '0', '0.2')
    sumi = 0
    for i in range(1,161):
    	previous_value = transfered_cp_x
    	current_value = current_cp[0]
    	sumi = sumi + (float(current_value[i])-float(previous_value[i]))**2
    l2array.append(sumi**0.5)

fig1 = plt.figure(1)
ax = plt.subplot(111)
ax.plot(modified_x_y[0],modified_x_y[1], label = 'Designed Airfoil')
ax.plot(target_airfoil[0],target_airfoil[1], label = 'Target Airfoil')
ax.plot(initial_airfoil[0],initial_airfoil[1], label = 'Initial Airfoil')
ax.legend()
fig1.savefig('result.png')                                                 #Saving Airfoils Plot

plot_target_cp = Xfoil('n0009', '5000000' , '0', '0.2')
plot_initial_cp = Xfoil('n0015', '5000000' , '0', '0.2')
plot_final_cp = Xfoil('Modified', '5000000' , '0', '0.2')

fig2 = plt.figure(2)
ax = plt.subplot(111)
ax.plot(plot_target_cp[1],plot_target_cp[0], label = 'Target Airfoil')
ax.plot(plot_initial_cp[1],plot_initial_cp[0], label = 'Initial Airfoil')
ax.plot(plot_final_cp[1],plot_final_cp[0], label = 'Designed Airfoil')
ax.legend()
fig2.savefig('cpresult.png')                                               #Saving Cp Plot

normx = []
for i in range(1,170):
	normx.append(i)

print(l2array)

fig2 = plt.figure(3)
ax = plt.subplot(111)
ax.plot(normx,l2array,'ro',label = 'L2 Norm')
ax.legend()
fig2.savefig('l2norm.png')  

