
# compare analytical solution to gekko solution

import numpy as np
import matplotlib.pyplot as plt
from gekko import*
from mpl_toolkits.mplot3d.axes3d import Axes3D

#analytical solution

def phi(x):
    phi = np.cos(x)
    return phi

def psi(x):
    psi = np.sin(2*x)
    return psi

def ua(x,t):
    #u = np.cos(x)*np.cos(3*t) + 1/6*np.sin(2*x)*np.sin(6*t)
    a = 18996.06 # ft/s speed of sound in steel
    c = a # 3 (from example problem)
    #u = 1/2*(np.cos(x-c*t)+np.cos(x+c*t)) - 1/(4*c)*(np.cos(2*(x+c*t)) -np.cos(2*(x-c*t))) 
    u = np.cos(x)*np.cos(a*t) + 1/(2*a)*np.sin(2*x)*np.sin(2*a*t)
    return u


# define time
tf = .0005
npt = 100#101
xf = 2*np.pi
npx = 100#40

time = np.linspace(0,tf,npt)
xpos = np.linspace(0,xf,npx)

for i in range(npx):
    usol = ua(xpos[i],time)
    if i == 0:
        ustora = usol
    else:
        ustora = np.vstack([ustora,usol])

for i in range(npt):
    if i ==0:
        xstor = xpos
    else:
        xstor = np.vstack([xstor,xpos])

for i in range(npx):
    if i == 0:
        tstor = time
    else:
        tstor = np.vstack([tstor,time])
        
xstor = xstor.T
    


#%%
# create gekko model 

m = GEKKO() # (remote=False) for local solution

m.time = time

x0 = phi(xpos) 
v0 = psi(xpos)
dx = xpos[1]-xpos[0]
a = 18996.06 # ft/s speed of sound in steel
c = m.Const(value = a)
dx = m.Const(value = dx)
u = [m.Var(value = x0[i]) for i in range(npx)]
v = [m.Var(value = v0[i]) for i in range(npx)]
#j = [m.Var() for i in range(npx)]

[m.Equation(u[i].dt()==v[i]) for i in range(npx)] 
# top difference eqution (forward) first order
#m.Equation(v[0].dt()==c**2*(1/dx**2)*(u[2]-2*u[1]+u[0]))
# second order
#m.Equation(v[0].dt()==c**2*(1/dx**2)*(-u[3] + 4*u[2])-5*u[1] +2*u[0])
# 
m.Equation(v[0].dt()==c**2 * (u[1] - 2.0*u[0] + u[npx-1])/dx**2 )
# central difference (middle)
[m.Equation(v[i+1].dt()== c**2 * (u[i+2] - 2.0*u[i+1] + u[i])/dx**2) for i in range(npx-2) ]
# bottom  (backward) first order
#m.Equation(v[npx-1].dt()==c**2*(1/dx**2)*(u[npx-1] -2*u[npx-2]+u[npx-3]))
# second order
m.Equation(v[npx-1].dt()== c**2 * (u[npx-2] - 2.0*u[npx-1] + u[0])/dx**2 )
# set options
#m.Equaiton(j == v[npx].dt())
m.options.imode = 4
m.options.solver = 1
m.options.nodes = 3

m.solve()
    
for i in range(npx):
    if i ==0:
        ustor = np.array([u[i]])
        tstor = np.array([m.time])
    else:
        ustor = np.vstack([ustor,u[i]])
        tstor = np.vstack([tstor,m.time])

for i in range(npt):
    if i == 0:
        xstor = xpos
    else:
        xstor = np.vstack([xstor,xpos])
xstor = xstor.T
t = tstor
ustor = np.array(ustor) 

#%%
# compute error
error = ustora - ustor

#%%
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1,projection='3d')
ax.set_xlabel('Distance (ft)', fontsize = 12)
ax.set_ylabel('Time (seconds)', fontsize = 12)
ax.set_zlabel('Position (ft)', fontsize = 12)
#plt.title('Analytical Solution')
ax.set_zlim((-1,1))
p = ax.plot_wireframe(xstor,tstor,ustora,rstride=1,cstride=1)
fig.savefig('analytical_3d.eps', dpi = 1200, Transparent = True)
fig.show()  

#%%
### WE USE THIS PLOT FOR THE GEKKO SOLUTION
# 3d gekko solution plot
fig = plt.figure(2)
ax = fig.add_subplot(1,1,1,projection='3d')
ax.set_xlabel('Distance (ft)', fontsize = 12)
ax.set_ylabel('Time (seconds)', fontsize = 12)
ax.set_zlabel('Position (ft)', fontsize = 12)
ax.set_zlim((-1,1))
#plt.title('GEKKO Solution')
p = ax.plot_wireframe(xstor,tstor,ustor,rstride=1,cstride=1)
fig.savefig('gekko_3d.eps', dpi = 1200, Transparent = True)
fig.show() 

#%%
### WE USE THIS PLOT FOR THE ANALYTICAL SOLUTION ###
# PLot analytical solution contour
plt.figure()  # start a new figure
plt.contour(xstor, tstor, ustora, 150)  # using 50 contour lines.
plt.colorbar()  # add a colorbar
plt.xlabel('X')  # labels for axes
plt.ylabel('Time')
plt.title('Analytical Solution')
plt.show()  # show plot

#%%
# plot gekko solution contour plot
plt.figure()  # start a new figure
plt.contour(xstor, tstor, ustor, 150)  # using 50 contour lines.
plt.colorbar()  # add a colorbar
plt.xlabel('X')  # labels for axes
plt.ylabel('Time')
plt.title('GEKKO Solution')
plt.show()  # show plot

#%%


# --- setup grid ---
#nx = 200  # number of points in x-direction
#ny = 150  # number of points in y-direction
#x = np.linspace(-5, 5, nx)  # nx points equally spaced between -5...5
#y = np.linspace(-6, 6, ny)  # ny points equally spaced between -6...6
#X, Y = np.meshgrid(x, y, indexing='ij')  # 2D array (matrix) of points across x and y
#Z = np.zeros((nx, ny))  # initialize output of size (nx, ny)

### WE USE THIS PLOT FOR THE DIFFERENCE BETWEEN THE ANALYTICAL AND GEKKO SOLUTIONS ###
# --- contour plot ---
plt.figure()  # start a new figure
plt.contour(xstor, tstor, error, 150)  # using 50 contour lines.
cbar = plt.colorbar() # add a colorbar
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Difference', fontsize=12)
plt.xlabel('Distance (ft)', fontsize = 12)  # labels for axes
plt.ylabel('Time (seconds)', fontsize = 12)
#plt.zlabel('Error')
#plt.title('Error (ft)')
plt.savefig('difference_gekko.eps', dpi = 1200, Transparent = True)
plt.show()  # show plot
 
#%%
per_error = error / ustora
per_error = np.abs(per_error)

# --- contour plot ---
plt.figure()  # start a new figure
plt.contour(xstor, tstor, per_error, 150)  # using 50 contour lines.
plt.colorbar(label = 'Error (%)')  # add a colorbar
plt.xlabel('X (ft)')  # labels for axes
plt.ylabel('Time (seconds)')
#plt.zlabel('Error (%)')
#plt.title('Error (ft)')
plt.savefig('error_gekko.eps')
plt.show()  # show plot







