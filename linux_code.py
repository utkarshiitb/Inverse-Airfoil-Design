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
from decimal import Decimal


run = "xfoil"

def Xfoil(name, Re, Alpha , Mach ):
    def Cmd(cmd):
        ps.stdin.write(cmd+'\n')
    try:
        os.remove(name+'.log')
    except :
        pass
    #    print ("no such file")
    # run xfoil
    ps = sp.Popen(run ,stdin=sp.PIPE,stderr=sp.PIPE,stdout=sp.PIPE, shell = True)
    ps.stderr.close()
    # comand part
    Cmd('load '+name+'.csv')
    Cmd('OPER')
    Cmd('visc '+str(Re))
    Cmd('mach'+str(Mach))
    Cmd('alfa' + str(Alpha))
    Cmd('cpwr')
    Cmd(name+'CP'+'.csv')
    Cmd(' ')     # escape OPER
    Cmd('quit')  # exit
    #resp = ps.stdout.read()
    #print "resp:",resp   # console ouput for debug
    ps.stdout.close()
    ps.stdin.close()
    ps.wait()
    #while (ps.returncode() == None):
    #    time.sleep(1)
    #ps.kill()
        
    #Getting X and CP in Array

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
    #ascii.write([xvalue, cpvalue], 'Modified.csv', names=['NACA6409', 'y'])
    #print cp_x
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
    #ascii.write([xvalues, yvalues], 'Modified.csv', names=['NACA0009', 'y'])
    #print xvalues
    #print yvalues

def create_xy(x_y_array):
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
            virtual_y.append(-abs((1/(slope[i]**2+1)*diff))**0.5)


#calculating the relaxation factor
    omega = []
    for i in range(len(target_cp)-1):
        omega.append(ep*(1/speed_sound)*(abs((target_press[i] - current_press[i])/rho
)**0.5))

#calculating displacement parallel to the y axis
    dy = []
    for i in range(len(target_cp)-1):
        dy.append(-omega[i]*virtual_y[i]*dt)
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


#airfoil_x = [1.0, 0.99249, 0.98146, 0.9685, 0.95462, 0.94039, 0.92603, 0.91162, 0.89718, 0.88271, 0.86823, 0.85372, 0.83919, 0.82464, 0.81007, 0.79548, 0.78088, 0.76626, 0.75162, 0.73697, 0.7223, 0.70762, 0.69293, 0.67822, 0.6635, 0.64878, 0.63404, 0.61929, 0.60454, 0.58978, 0.57501, 0.56023, 0.54546, 0.53068, 0.51589, 0.50111, 0.48633, 0.47154, 0.45676, 0.44199, 0.42723, 0.41251, 0.39789, 0.38349, 0.36915, 0.35486, 0.34058, 0.32633, 0.3121, 0.2979, 0.28374, 0.26961, 0.25553, 0.2415, 0.22752, 0.2136, 0.19975, 0.18598, 0.1723, 0.15871, 0.14524, 0.1319, 0.1187, 0.10569, 0.0929, 0.08038, 0.06823, 0.0566, 0.04573, 0.03597, 0.02767, 0.02098, 0.01575, 0.0117, 0.00856, 0.00609, 0.00416, 0.00265, 0.0015, 0.00069, 0.00019, 0.0, 0.00012, 0.00061, 0.00152, 0.00285, 0.00458, 0.0067, 0.00924, 0.01228, 0.01595, 0.02048, 0.02614, 0.03333, 0.04237, 0.05333, 0.06587, 0.07946, 0.09369, 0.10831, 0.12318, 0.13822, 0.15339, 0.16859, 0.18376, 0.19887, 0.21394, 0.22896, 0.24395, 0.2589, 0.27383, 0.28874, 0.30363, 0.3185, 0.33336, 0.3482, 0.36304, 0.37788, 0.39274, 0.40777, 0.42296, 0.43817, 0.45339, 0.46859, 0.4838, 0.49899, 0.51417, 0.52935, 0.54452, 0.55968, 0.57484, 0.58999, 0.60514, 0.62028, 0.63541, 0.65055, 0.66567, 0.6808, 0.69592, 0.71103, 0.72615, 0.74126, 0.75636, 0.77147, 0.78657, 0.80167, 0.81677, 0.83186, 0.84696, 0.86205, 0.87714, 0.89223, 0.90731, 0.92239, 0.93744, 0.95238, 0.96698, 0.9806, 0.99216, 1.0]

#airfoil_y=[9.45E-04,3.22E-03,6.53E-03,1.03E-02,1.43E-02,1.83E-02,2.23E-02,2.61E-02,2.99E-02,3.36E-02,3.72E-02,4.08E-02,4.42E-02,4.75E-02,5.08E-02,5.40E-02,5.71E-02,6.01E-02,6.30E-02,6.58E-02,6.85E-02,7.11E-02,7.37E-02,7.61E-02,7.85E-02,8.08E-02,8.29E-02,8.50E-02,8.70E-02,8.88E-02,9.06E-02,9.23E-02,9.39E-02,9.53E-02,9.67E-02,9.79E-02,9.91E-02,0.1001306,0.1010541,0.1018632,0.1025559,0.1031294,0.103578,0.1038499,0.1039143,0.1037726,0.1034228,0.1028628,0.1020904,0.1011033,9.99E-02,9.85E-02,9.68E-02,9.50E-02,9.29E-02,9.05E-02,8.80E-02,8.52E-02,8.22E-02,7.89E-02,7.54E-02,7.16E-02,6.76E-02,6.33E-02,5.88E-02,5.41E-02,4.91E-02,4.39E-02,3.86E-02,3.34E-02,2.86E-02,2.42E-02,2.05E-02,1.72E-02,1.44E-02,1.19E-02,9.62E-03,7.52E-03,5.55E-03,3.68E-03,1.90E-03,0.00E00,-1.42E-03,-3.08E-03,-4.67E-03,-6.12E-03,-7.41E-03,-8.55E-03,-9.56E-03,-1.05E-02,-1.12E-02,-1.19E-02,-1.24E-02,-1.27E-02,-1.28E-02,-1.25E-02,-1.18E-02,-1.07E-02,-9.48E-03,-8.04E-03,-6.47E-03,-4.83E-03,-3.15E-03,-1.46E-03,2.12E-04,1.85E-03,3.44E-03,4.97E-03,6.43E-03,7.82E-03,9.14E-03,1.04E-02,1.15E-02,1.26E-02,1.35E-02,1.44E-02,1.51E-02,1.57E-02,1.63E-02,1.67E-02,1.71E-02,1.75E-02,1.78E-02,1.81E-02,1.84E-02,1.86E-02,1.88E-02,1.90E-02,1.91E-02,1.91E-02,1.92E-02,1.91E-02,1.91E-02,1.90E-02,1.88E-02,1.86E-02,1.83E-02,1.80E-02,1.76E-02,1.72E-02,1.67E-02,1.62E-02,1.56E-02,1.50E-02,1.43E-02,1.36E-02,1.28E-02,1.19E-02,1.10E-02,1.01E-02,9.09E-03,8.04E-03,6.93E-03,5.78E-03,4.58E-03,3.33E-03,2.07E-03,8.55E-04,-2.09E-04, 9.45E-04]

#target_cp = [0.18694, 0.17453418204182033, 0.15485056603773595, 0.13073023112480753, 0.10425132780082987, 0.07694614239482195, 0.04954621013133212, 0.0224777900552486, -0.004088472306756, -0.030084470018170873, -0.05545767835550194, -0.08026905853952924, -0.10452759348612795, -0.1282305790108564, -0.15152683544303805, -0.17449556359252572, -0.19718363636363637, -0.2196808313253011, -0.24205333534015658, -0.2643431125827815, -0.2866376339554485, -0.3089802227573751, -0.33143765803732683, -0.3540753461770017, -0.37694638554216875, -0.40007546987951803, -0.42354354430379737, -0.4473887823990355, -0.47163920337960175, -0.49637973445986733, -0.5216946409173205, -0.5476533534743202, -0.574233107617896, -0.6015364307320024, -0.6296289824348881, -0.658540952092177, -0.6883628901032179, -0.7191928267477202, -0.751074692635423, -0.7840866992068334, -0.8185507138499085, -0.8543201832620648, -0.8913752264381885, -0.9295235460122699, -0.969289249692497, -1.0108576942046856, -1.0545500989486705, -1.1005455707196028, -1.1492206451612903, -1.2013217061021169, -1.2569161038148842, -1.316634418604651, -1.381127364096081, -1.4514401593371575, -1.5290287057308434, -1.6162152941176469, -1.7173132083054252, -1.8384765217391303, -1.9743324223602483, -2.121246638115632, -2.2758354996299035, -2.4349664620107445, -2.5957983, -2.7531114759825326, -2.9010176405529955, -3.036961705882353, -3.202760221518987, -3.3839000690448793, -3.4896651851851854, -3.6029318387096776, -3.720716976744186, -3.8378843820224717, -3.9576501612903225, -4.071518154506438, -4.178956237623763, -4.269960114942529, -4.335152905405405, -4.35145297029703, -4.29477, -4.1367112068965515, -3.8473608333333336, -3.34402, -2.8642955172413793, -2.0863607142857146, -1.2200202298850573, -0.4186895689655171, 0.16957005714285717, 0.5576544827586207, 0.7996776086956522, 0.9385185769230769, 1.006982218430034, 1.0219082085561497, 1.0080328571428572, 0.9739767417677643, 0.92791, 0.8769760512820512, 0.8278147878787879, 0.7815886776859504, 0.7386159756097561, 0.6980901339285714, 0.6591407063197027, 0.6211469003690037, 0.583737435105068, 0.5470137039350407, 0.5113710542929293, 0.4776680661577608, 0.4482087086513995, 0.42274629770992367, 0.3999376425855513, 0.37932972222222217, 0.3605913316582914, 0.34347729205753597, 0.32780413707165107, 0.3134192369727047, 0.30017519777503093, 0.287957120838471, 0.27665198402948404, 0.26618060698957696, 0.2564445476772616, 0.24731335570469795, 0.23874986577181206, 0.2308041047503045, 0.22339996354799516, 0.21649500000000002, 0.2100420411871593, 0.20401321234119782, 0.19836778247734138, 0.19307439951719976, 0.18810905364677516, 0.1834470481927711, 0.1790706016847172, 0.1749552134696332, 0.1710834455802766, 0.1674684314903846, 0.16407260504201682, 0.16088594594594594, 0.1579054349130174, 0.15511436112777446, 0.15252194244604317, 0.1501316366906475, 0.14793480215827337, 0.1459374220623501, 0.14415338129496402, 0.14260433453237412, 0.14130888422315538, 0.14028827338129496, 0.13956331332533012, 0.1391744537815126, 0.13916953753753752, 0.1396160950661853, 0.14059896925858953, 0.1422254479418886, 0.14465529016493583, 0.14812126865671643, 0.15295074626865673, 0.15963869285254342, 0.1682481911966987, 0.17737652873563217, 0.18419887530562346, 0.18694]


target_cp_x = Xfoil('n6409', '5000000' , '0', '0.3')
modified_x_y = ReadXY('Modified')
transfered_cp_x = get_cp_on_x(target_cp_x,modified_x_y[0])
for x in xrange(1,10000):
    plt.clf()
    current_cp = Xfoil('Modified', '5000000', '0', '0.3')
    modified_x_y = algorithm(transfered_cp_x, current_cp[0], modified_x_y, 0.3)
    create_xy(modified_x_y)
    #plt.plot(modified_x_y[0], modified_x_y[1])
    #plt.show()




#current_cp = cp_x[0]
#mach = 0.3
#airfoil_xy = [0, 0]
#airfoil_xy[0] = airfoil_x
#airfoil_xy[1] = airfoil_y
#letsc = get_cp_on_x(cp_x, airfoil_x)
#print letsc
#new_airfoil = raja_bene_wahi(target_cp, current_cp, airfoil_xy, mach)

#plt.plot(airfoil_x, airfoil_y)
#plt.plot(new_airfoil[0], new_airfoil[1])

#plt.show()

