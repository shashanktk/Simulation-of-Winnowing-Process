#!/usr/bin/env python
# coding: utf-8

# In[6]:


#Shashank Tumkur Karnick #matriculation = 229777

#Euler Method

import numpy as np
import math
from matplotlib import pyplot as plt

#Step Size
dt = 0.00001

#start time
tstart = 0

#End Time
tend = 20

#Number of iterations
N = int((tend-tstart)/dt)

#Time vector
t = np.linspace(tstart,tend,N)

#Given conditions(Dimesnions and Flow Properties)
x0 = 0.5
y0 = 0.5
xc = 0.55
yc = -0.5


# dimesions vectors for heavy particles and light particles

x_g    = np.zeros(len(t))     #x_g direction for grains
x_g[0] = x0
y_g    = np.zeros(len(t))     #y_g direction for grains
y_g[0] = y0
x_c    = np.zeros(len(t))     #x_c direction for chaffs
x_c[0] = x0
y_c    = np.zeros(len(t))     #y_c direction for chaffs
y_c[0] = y0


#Given conditions

#Fluid properties
h     = 0.1
u0    = 0.2
rho_f = 1.2                     #Density of air
Mue   = 1.8*10**-5              #Viscocity of air
 
#Grain particle properties
rho_g = 750                     #Density of grain particle
dia_g = 0.0025                  #Diameter of grain particles

#Chaff particle properties
rho_c = 50                      #density of chaff prticle
dia_c = 0.00325                 #diameter of chaff particles

#Air velocity function
ufg_x = np.zeros(len(t))
ufg_y = np.zeros(len(t))
ufc_x = np.zeros(len(t))
ufc_y = np.zeros(len(t))

#Reynold number and Cd vectors
Re_g  = np.zeros(len(t))
Re_c  = np.zeros(len(t))
Cd_g  = np.zeros(len(t))
Cd_c  = np.zeros(len(t))

#Drag force vector
Fdrag_g_x = np.zeros(len(t))
Fdrag_c_x = np.zeros(len(t))
Fdrag_g_y = np.zeros(len(t))
Fdrag_c_y = np.zeros(len(t))

#Particle velocity function
upg_x = np.zeros(len(t))          # Grain velocity function in x direction
upg_y = np.zeros(len(t))          # Grain velocity function in y direction
upc_x = np.zeros(len(t))          # Chaff velocity function in x direction
upc_y = np.zeros(len(t))          # Chaff velocity function in y direction

#Volume of grain particle and chaff particle
Vol_g = (math.pi*dia_g**3)/6          # Grain particle
Vol_c = (math.pi*dia_c**3)/6          # chaff particle

#Forces involved
F_gravity_g_y  = -rho_g*Vol_g*9.8    #Gravitational force for grain particles
F_gravity_c_y  = -rho_c*Vol_c*9.8    #Gravitational force for chaff particles
F_gravity_g_x  = 0                   #Gravitational force for grain particles
F_gravity_c_x  = 0                   #Gravitational force for grain particles
F_buoyancy_g_y = rho_f*Vol_g*9.8     #Gravitational force for grain particles
F_buoyancy_c_y = rho_f*Vol_c*9.8     #Gravitational force for chaff particles
F_buoyancy_g_x = 0                   #Gravitational force for grain particles
F_buoyancy_c_x = 0                   #Gravitational force for grain particles

#Function handles and Computational loops for Euler Method

#For x_g

def f_g_x(t,y):
    return y

def f_up_g_x(t,y):
    return y/(rho_g*Vol_g)



# For y_g

def f_g_y(t,y):
    return y
def f_upg_y(t,y):
    return F_gravity_g_y/(rho_g*Vol_g) + F_buoyancy_g_y/(rho_g*Vol_g)+ y/(rho_g*Vol_g)


#For x_c

def f_c_x(t,y):
    return y
def f_upc_x(t,y):
    return y/(rho_c*Vol_c)


#For y_c

def f_c_y(t,y):
    return y
def f_upc_y(t,y):
    return F_gravity_c_y/(rho_c*Vol_c) + F_buoyancy_c_y/(rho_c*Vol_c) + y/(rho_c*Vol_c)


#Loop

for i in range(0,N-1):
    ufg_x[i] = 6.2*u0*math.sqrt((h/x_g[i]))*math.exp(-50*(y_g[i])**2/(x_g[i])**2)
    ufc_x[i] = 6.2*u0*math.sqrt((h/x_c[i]))*math.exp(-50*(y_c[i])**2/(x_c[i])**2)
    ufg_y[i] = 0
    ufc_y[i] = 0
    Re_g[i]  = rho_f*math.sqrt((ufg_x[i]-upg_x[i])**2+(ufg_y[i]-upg_y[i])**2)*dia_g/Mue
    Re_c[i]  = rho_f*math.sqrt((ufc_x[i]-upc_x[i])**2+(ufc_y[i]-upc_y[i])**2)*dia_c/Mue

    if Re_g[i] < 800:
        Cd_g[i] = (24/Re_g[i])*(1+0.15*Re_g[i]**0.687)
        Cd_c[i] = (24/Re_c[i])*(1+0.15*Re_c[i]**0.687)
    else:
        Cd_g[i] = 0.44
        Cd_c[i] = 0.44

    Fdrag_g_x[i] = math.pi*0.5*dia_g**2*rho_f*Cd_g[i]*math.sqrt((ufg_x[i]-upg_x[i])**2+(ufg_y[i]-upg_y[i])**2)*(ufg_x[i]-upg_x[i])
    Fdrag_c_x[i] = math.pi*0.5*dia_c**2*rho_f*Cd_c[i]*math.sqrt((ufc_x[i]-upc_x[i])**2+(ufc_y[i]-upc_y[i])**2)*(ufc_x[i]-upc_x[i])
    Fdrag_g_y[i] = math.pi*0.5*dia_g**2*rho_f*Cd_g[i]*math.sqrt((ufg_x[i]-upg_x[i])**2+(ufg_y[i]-upg_y[i])**2)*(ufg_y[i]-upg_y[i])
    Fdrag_c_y[i] = math.pi*0.5*dia_c**2*rho_f*Cd_c[i]*math.sqrt((ufc_x[i]-upc_x[i])**2+(ufc_y[i]-upc_y[i])**2)*(ufc_y[i]-upc_y[i])

        #computational loop for grain(x-axis)
    x_g[i+1]   = x_g[i]   + dt *f_g_x(t[i],upg_x[i])
    upg_x[i+1] = upg_x[i] + dt *f_up_g_x(t[i],Fdrag_g_x[i])
        
        #computational loop for grain(y-axis)
    y_g[i+1]   = y_g[i]   + dt *f_g_y(t[i],upg_y[i])
    upg_y[i+1] = upg_y[i] + dt *f_upg_y(t[i],Fdrag_g_y[i])
        
        #computational loop for chaff (x-axis)
    x_c[i+1]   = x_c[i]   + dt *f_c_x(t[i],upc_x[i])
    upc_x[i+1] = upc_x[i] + dt *f_upc_x(t[i],Fdrag_c_x[i])

        #computational loop for chaff(y-axis)
    y_c[i+1]   = y_c[i]   + dt *f_c_y(t[i],upc_y[i])
    upc_y[i+1] = upc_y[i] + dt *f_upc_y(t[i],Fdrag_c_y[i])
        




#Plots
plt.xlim(0.4,0.7)
plt.ylim(-0.5,0.6)
plt.plot(x_g,y_g,'green',label='Grain')
plt.plot(x_c,y_c,'blue',label='Chaff')
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.legend()
plt.show()


# In[8]:


#Runge-kutta method


import numpy as np
import math
from matplotlib import pyplot as plt

#Step Size
dt = 0.00001

#start time
tstart = 0

#End Time
tend = 20

#Number of iterations
N = int((tend-tstart)/dt)

#Time vector
t = np.linspace(tstart,tend,N)

#Given conditions(Dimesnions and Flow Properties)
x0 = 0.5
y0 = 0.5
xc = 0.55
yc = -0.5


# dimesions vectors for heavy particles and light particles

x_g    = np.zeros(len(t))     #x_g direction for grains
x_g[0] = x0
y_g    = np.zeros(len(t))     #y_g direction for grains
y_g[0] = y0
x_c    = np.zeros(len(t))     #x_c direction for chaffs
x_c[0] = x0
y_c    = np.zeros(len(t))     #y_c direction for chaffs
y_c[0] = y0


#Given conditions

#Fluid properties
h     = 0.1
u0    = 0.2
rho_f = 1.2                     #Density of air
Mue   = 1.8*10**-5              #Viscocity of air
 
#Grain particle properties
rho_g = 750                     #Density of grain particle
dia_g = 0.0025                  #Diameter of grain particles

#Chaff particle properties
rho_c = 50                      #density of chaff prticle
dia_c = 0.00325                 #diameter of chaff particles

#Air velocity function
ufg_x = np.zeros(len(t))
ufg_y = np.zeros(len(t))
ufc_x = np.zeros(len(t))
ufc_y = np.zeros(len(t))

#Reynold number and Cd vectors
Re_g  = np.zeros(len(t))
Re_c  = np.zeros(len(t))
Cd_g  = np.zeros(len(t))
Cd_c  = np.zeros(len(t))

#Drag force vector
Fdrag_g_x = np.zeros(len(t))
Fdrag_c_x = np.zeros(len(t))
Fdrag_g_y = np.zeros(len(t))
Fdrag_c_y = np.zeros(len(t))

#Particle velocity function
upg_x = np.zeros(len(t))          # Grain velocity function in x direction
upg_y = np.zeros(len(t))          # Grain velocity function in y direction
upc_x = np.zeros(len(t))          # Chaff velocity function in x direction
upc_y = np.zeros(len(t))          # Chaff velocity function in y direction

#Volume of grain particle and chaff particle
Vol_g = (math.pi*dia_g**3)/6          # Grain particle
Vol_c = (math.pi*dia_c**3)/6          # chaff particle

#Forces involved
F_gravity_g_y  = -rho_g*Vol_g*9.8    #Gravitational force for grain particles
F_gravity_c_y  = -rho_c*Vol_c*9.8    #Gravitational force for chaff particles
F_gravity_g_x  = 0                   #Gravitational force for grain particles
F_gravity_c_x  = 0                   #Gravitational force for grain particles
F_buoyancy_g_y = rho_f*Vol_g*9.8     #Gravitational force for grain particles
F_buoyancy_c_y = rho_f*Vol_c*9.8     #Gravitational force for chaff particles
F_buoyancy_g_x = 0                   #Gravitational force for grain particles
F_buoyancy_c_x = 0                   #Gravitational force for grain particles

#Function handles and Computational loops for Euler Method

#For x_g

def f_g_x(t,y):
    return y

def f_up_g_x(t,y):
    return y/(rho_g*Vol_g)



# For y_g

def f_g_y(t,y):
    return y
def f_upg_y(t,y):
    return F_gravity_g_y/(rho_g*Vol_g) + F_buoyancy_g_y/(rho_g*Vol_g)+ y/(rho_g*Vol_g)


#For x_c

def f_c_x(t,y):
    return y
def f_upc_x(t,y):
    return y/(rho_c*Vol_c)


#For y_c

def f_c_y(t,y):
    return y
def f_upc_y(t,y):
    return F_gravity_c_y/(rho_c*Vol_c) + F_buoyancy_c_y/(rho_c*Vol_c) + y/(rho_c*Vol_c)


#Loop

for i in range(0,N-1):
    ufg_x[i] = 6.2*u0*math.sqrt((h/x_g[i]))*math.exp(-50*(y_g[i])**2/(x_g[i])**2)
    ufc_x[i] = 6.2*u0*math.sqrt((h/x_c[i]))*math.exp(-50*(y_c[i])**2/(x_c[i])**2)
    ufg_y[i] = 0
    ufc_y[i] = 0
    Re_g[i]  = rho_f*math.sqrt((ufg_x[i]-upg_x[i])**2+(ufg_y[i]-upg_y[i])**2)*dia_g/Mue
    Re_c[i]  = rho_f*math.sqrt((ufc_x[i]-upc_x[i])**2+(ufc_y[i]-upc_y[i])**2)*dia_c/Mue

    if Re_g[i] < 800:
        Cd_g[i] = (24/Re_g[i])*(1+0.15*Re_g[i]**0.687)
        Cd_c[i] = (24/Re_c[i])*(1+0.15*Re_c[i]**0.687)
    else:
        Cd_g[i] = 0.44
        Cd_c[i] = 0.44

    Fdrag_g_x[i] = math.pi*0.5*dia_g**2*rho_f*Cd_g[i]*math.sqrt((ufg_x[i]-upg_x[i])**2+(ufg_y[i]-upg_y[i])**2)*(ufg_x[i]-upg_x[i])
    Fdrag_c_x[i] = math.pi*0.5*dia_c**2*rho_f*Cd_c[i]*math.sqrt((ufc_x[i]-upc_x[i])**2+(ufc_y[i]-upc_y[i])**2)*(ufc_x[i]-upc_x[i])
    Fdrag_g_y[i] = math.pi*0.5*dia_g**2*rho_f*Cd_g[i]*math.sqrt((ufg_x[i]-upg_x[i])**2+(ufg_y[i]-upg_y[i])**2)*(ufg_y[i]-upg_y[i])
    Fdrag_c_y[i] = math.pi*0.5*dia_c**2*rho_f*Cd_c[i]*math.sqrt((ufc_x[i]-upc_x[i])**2+(ufc_y[i]-upc_y[i])**2)*(ufc_y[i]-upc_y[i])

        #computational loop for grain(x-axis)
    k1= dt* f_g_x(t[i],upg_x[i])
    k2= dt* f_g_x(t[i],upg_x[i]+k1/2)
    k3= dt* f_g_x(t[i],upg_x[i]+k2/2)
    k4= dt* f_g_x(t[i],upg_x[i]+k3)
    x_g[i+1]   = x_g[i]   + (k1+2*k2+2*k3+k4)/6


    k1= dt* f_up_g_x(t[i],Fdrag_g_x[i])
    k2= dt* f_up_g_x(t[i],Fdrag_g_x[i]+k1/2)
    k3= dt* f_up_g_x(t[i],Fdrag_g_x[i]+k2/2)
    k4= dt* f_up_g_x(t[i],Fdrag_g_x[i]+k3)
    upg_x[i+1] = upg_x[i] + (k1+2*k2+2*k3+k4)/6
        
        #computational loop for grain(y-axis)

    k1= dt* f_g_y(t[i],upg_y[i])
    k2= dt* f_g_y(t[i],upg_y[i]+k1/2)
    k3= dt* f_g_y(t[i],upg_y[i]+k2/2)
    k4= dt* f_g_y(t[i],upg_y[i]+k3)
    y_g[i+1]   = y_g[i]   + (k1+2*k2+2*k3+k4)/6


    k1= dt* f_upg_y(t[i],Fdrag_g_y[i])
    k2= dt* f_upg_y(t[i],Fdrag_g_y[i]+k1/2)
    k3= dt* f_upg_y(t[i],Fdrag_g_y[i]+k2/2)
    k4= dt* f_upg_y(t[i],Fdrag_g_y[i]+k3)
    upg_y[i+1] = upg_y[i] + (k1+2*k2+2*k3+k4)/6
    
        
        #computational loop for chaff (x-axis)
    k1= dt* f_c_x(t[i],upc_x[i])
    k2= dt* f_c_x(t[i],upc_x[i]+k1/2)
    k3= dt* f_c_x(t[i],upc_x[i]+k2/2)
    k4= dt* f_c_x(t[i],upc_x[i]+k3)
    x_c[i+1]   = x_c[i]   + (k1+2*k2+2*k3+k4)/6

    k1= dt* f_upc_x(t[i],Fdrag_c_x[i])
    k2= dt* f_upc_x(t[i],Fdrag_c_x[i]+k1/2)
    k3= dt* f_upc_x(t[i],Fdrag_c_x[i]+k2/2)
    k4= dt* f_upc_x(t[i],Fdrag_c_x[i]+k3)
    upc_x[i+1]   = upc_x[i]   + (k1+2*k2+2*k3+k4)/6


        #computational loop for chaff(y-axis)
    k1= dt* f_c_y(t[i],upc_y[i])
    k2= dt* f_c_y(t[i],upc_y[i]+k1/2)
    k3= dt* f_c_y(t[i],upc_y[i]+k2/2)
    k4= dt* f_c_y(t[i],upc_y[i]+k3)
    y_c[i+1]   = y_c[i]   + (k1+2*k2+2*k3+k4)/6

    k1= dt* f_upc_y(t[i],Fdrag_c_y[i])
    k2= dt* f_upc_y(t[i],Fdrag_c_y[i]+k1/2)
    k3= dt* f_upc_y(t[i],Fdrag_c_y[i]+k2/2)
    k4= dt* f_upc_y(t[i],Fdrag_c_y[i]+k3) 

    upc_y[i+1] = upc_y[i] + (k1+2*k2+2*k3+k4)/6
        




#Plots
plt.xlim(0.4,0.7)
plt.ylim(-0.5,0.6)
plt.plot(x_g,y_g,'green',label='grain')
plt.plot(x_c,y_c,'blue',label='chaff')
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.legend()
plt.show()


# In[ ]:




