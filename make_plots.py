# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 10:03:15 2021

@author: 00ery
"""
#LIBRARIES
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.patches as mpatches

#increasing pixel density (300 should give very good quality)
sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})

#writing contents of file to list (for floating point numbers)  
def file_to_list(file_name):
    f = open(file_name, 'r')
    f_lines = f.readlines()
    for i in range(len(f_lines)):
        f_lines[i] = float(f_lines[i])
    return f_lines

'''
x = file_to_list("nc1_bonded_fraction.txt")
x = x[1:len(x)]
x_smooth = []
for i in range(0,len(x),2):
    x_smooth.append(x[i]/2+x[i+1]/2)
plt.plot(x_smooth[59:])
'''

#PLOTTING SOMETHING WITH ERROR BARS

'''
cf = []
for i in range(100,200,10):
	cf.append(i/100)
conc = []
for i in cf:
    conc.append(1/(4*4*16*(i**3)))
'''   
    
    





    
y1 = file_to_list("OP.txt")
#y2 = file_to_list("broken100200.txt")


'''
#Savitzky-Golay Filter
from scipy.signal import savgol_filter
y1 = savgol_filter(y1,15,2)
y2 = savgol_filter(y2,15,2)
'''

x=np.linspace(60,95.7,72)


#converting to dataframe    
data1 = {"Concentration":x,
        "Quantity":y1}
data1 = pd.DataFrame(data1)

'''
data2 = {"Concentration":x,
        "Quantity":y2}
data2 = pd.DataFrame(data2)
'''






#making plot
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_style("darkgrid")
plt.ticklabel_format(style="plain")
sns.lineplot(data=data1, x="Concentration", y="Quantity", color="#8A2BE2")
#sns.lineplot(data=data2, x="Concentration", y="Quantity", color="#00FF00",label='broken')




plt.xlabel("Time $(10^6 \u03C4)$", fontsize=15)
plt.ylabel(r"$\langle  cos^2 \theta   \rangle$", fontsize=15)



plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

   
#error bars
error1 = np.array(file_to_list("OP_STDEV.txt"))/(1024**0.5)
plt.fill_between(x, y1-error1, y1+error1, alpha=0.15, facecolor='#8A2BE2')
                 


#DECORATIONS

plt.axvline(60, c='r', linestyle='dashed')
plt.axvline(65.7, c='r', linestyle='dashed')
plt.text(60.7,0.675,'Applying \nstrain \n(0 to 3)', fontsize="15", c='r')

#plt.plot(66.1,0.72,'ro', c='black') 
#plt.text(60.8,0.713,'yield point', fontsize="12", c='black')




'''
#MAKING NICE PLOTS FROM DATA FROM A FILE 

#converting to dataframe (plus some preprocessing)
y1 = np.array(file_to_list("force_CF1.8.txt"))
y2 = np.array(file_to_list("force_CF1.5.txt"))
y3 = np.array(file_to_list("force_CF1.4.txt"))
y4 = np.array(file_to_list("force_CF1.2.txt"))
y5 = np.array(file_to_list("force_CF1.0.txt"))

#FOR STRESS STRAIN
#making all numbers positive and subtracting baseline
y1 = [abs(x) for x in y1]
y2 = [abs(x) for x in y2]
y3 = [abs(x) for x in y3]
y4 = [abs(x) for x in y4]
y5 = [abs(x) for x in y5]

baseline = (np.mean(y1[3:240])+np.mean(y2[3:240])+np.mean(y3[3:240])+np.mean(y4[3:240])+np.mean(y5[3:240]))/5

y1=y1[241:] - baseline
y2=y2[241:] - baseline
y3=y3[241:] - baseline
y4=y4[241:] - baseline
y5=y5[241:] - baseline
'''

'''
#interpolating to de-noise
y1_new=[]
y2_new=[]
for i in range(0,160,2):
    y1_new.append((y1[i]+y1[i+1])/2)
    y2_new.append((y2[i]+y2[i+1])/2)
y1 = y1_new
y2 = y2_new
'''
'''
#calculating actual stress i.e. force/area
z = np.linspace(64,448,40)
area = 262144/z
y1 = y1/(area*1.8*1.8)
y2 = y2/(area*1.5*1.5)
y3 = y3/(area*1.4*1.4)
y4 = y4/(area*1.2*1.2)
y5 = y5/(area*1*1)



data1 = {"Timestep":np.linspace(0,6,40),
        "Quantity":y1}
data1 = pd.DataFrame(data1)

data2 = {"Timestep":np.linspace(0,6,40),
        "Quantity":y2}
data2 = pd.DataFrame(data2)

data3 = {"Timestep":np.linspace(0,6,40),
        "Quantity":y3}
data3 = pd.DataFrame(data3)

data4 = {"Timestep":np.linspace(0,6,40),
        "Quantity":y4}
data4 = pd.DataFrame(data4)

data5 = {"Timestep":np.linspace(0,6,40),
        "Quantity":y5}
data5 = pd.DataFrame(data5)


#making plot
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_style("darkgrid")
plt.ticklabel_format(style="plain")
sns.lineplot(data=data1, x="Timestep", y="Quantity", label="$W=1.3×10^7$", color="#0000FF", linewidth=2)
sns.lineplot(data=data2, x="Timestep", y="Quantity", label="$W=1.1×10^7$", color="#00FF00", linewidth=2)
sns.lineplot(data=data3, x="Timestep", y="Quantity", label="$W=1.0×10^7$", color="#FFFF00", linewidth=2)
sns.lineplot(data=data4, x="Timestep", y="Quantity", label="$W=8.2×10^6$", color="#FFA500", linewidth=2)
sns.lineplot(data=data5, x="Timestep", y="Quantity", label="$W=6.8×10^6$", color="#ff0000", linewidth=2)




                 
plt.xlabel("Strain", fontsize=15)
plt.ylabel("Stress ($\sigma \u03C4^{-4}	$)", fontsize=15)
#plt.title("Stress-strain curve ($\dot{\epsilon} = 4.375×10^{-5} \u03C4^{-1}$)", fontsize=20)

plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
        
#plt.legend()
'''


'''
#CALCULATING WORK DONE
stress_diff = y1-y2
work = np.sum(stress_diff)*(7/40)
print(work)

#Some decorations
plt.axhline(0, c='g', linestyle='dashed')
plt.text(6,-0.1,'No bonds broken', fontsize="15", c='g')
#plt.axvline(100, c='g', linestyle='dashed')
#plt.text(122,0.2,'Compressing box', fontsize="15", c='g')

rect=mpatches.Rectangle((100,-0.5),20,5, 
                        fill=False,
                        alpha=0.1,
                        color="purple",
                       linewidth=2,    
                       facecolor="green")
#plt.gca().add_patch(rect)
'''