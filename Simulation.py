# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 13:20:42 2018

@author: Brandon
"""

#from __future__ import division # compatibility with python 2.7
from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d.axes3d import Axes3D
#%%

# make an empty dictionary to store values
res = {}

# start time
start_time = time.time()

#BuildModle####################################################################
sim = GEKKO()


#Define Timespace
tf  = 1 # sec, length of simulation 
npt = 30 # number of time discretizations
nit = 1 # number of itterations to solve 
ti = np.linspace(0,tf,npt) # times for plotting
                        
sim.time = np.linspace(0,tf,npt)       

# Define Rod Discretizations
TVD_d = 4800 # ft, lenth of rod
npx = 30 # number of rod discretizations
dx = TVD_d/(npx-1) # #ft lenth of rod discretizations
xi = np.linspace(0,TVD_d,npx) # possitions allong rod (for plotting)

################################
# Conventional Rod Pump Unit Geometry
# API geometry dimension values
Ag=210.0
Cg=120.3
Ig=120.0
Pg=148.5
Hg=237.88
Gg=86.88
Rg=47.0

#lengths from FIG. 1 - Beam Pumping Unit Shown as a Four-Bar Linkage
L_1 = Rg
L_2 = np.sqrt((Hg-Gg)**2.0+Ig**2.0)
L_3 = Cg
L_4 = Pg
L_5 = Ag

#Simulation########################################

#Constants
sim.API       = sim.Const(value = 45)               #API gravity of fluid, unitless
sim.c         = sim.Const(value = 0.000013)           #Compressibility, psi^-1
sim.k         = sim.Const(value = 15)                 #Permeability, md
sim.Bo        = sim.Const(value = 1.2)               #FVF, rb/STB
sim.A_d       = sim.Const(value = 2)                #Drainage Area, Acres
sim.sw        = sim.Const(value = 0.2)               #Water Saturation
sim.porosity  = sim.Const(value = 0.08)        #Porosity, unitless
sim.gamma_E   = sim.Const(value = 1.78)         #Euler Constant
sim.C_a       = sim.Const(value = 31.6)             #Drainage Area Shape Factor (Circular)
sim.rw        = sim.Const(value = 0.328)             #Welbore radius, ft
sim.S         = sim.Const(value = 0)                  #unitless
sim.u_visc    = sim.Const(value = 1.5)                # Viscosity, cp
sim.h_pz      = sim.Const(value = 8)              #pay zone thickness, ft
sim.D_t       = sim.Const(value = 2.5)              # tubing diameter, in
sim.St_length = sim.Const(value = 85)         # rod pump stroke length, in

sim.g     = sim.Const(value = 32.2)             # acceleration due to gravity, ft/s^3
sim.g_conv= sim.Const(value = 32.2)           # lbf conversion , lb-ft/s^2-lbf
sim.rho_r = sim.Const(value = 490)              # lbs/ft^3, density of rod steel
sim.rho_w = sim.Const(value = 62.3 )          # lbs/ft^3, density of water at standard conditions
sim.a     = sim.Const(value =18996.06 )          # ft/s speed of sound in steel
sim.D_r   = sim.Const(value = 1.0)               # in, diameter of rod string
sim.Ac    = sim.Const(value= sim.D_r.value**2/4.0*np.pi) # in^2, cross sectional area of rod 
sim.nu    = sim.Const(value = 0.01)              # unitless, damping coefficient
sim.pi    = sim.Const(value=np.pi)
sim.E     = sim.Const(value = 32025000.0)        # psi sucker rod modulus of elasticity
sim.alpha = sim.Const(value = 0.0)            # pump parameter, unitless
sim.beta  = sim.Const(value = 1.0)            # pump parameter, unitless 

sim.L_1 = sim.Const(value =L_1) # unit geometry 
sim.L_2 = sim.Const(value =L_2)  # unit geometry 
sim.L_3 = sim.Const(value =L_3)  # unit geometry 
sim.L_4 = sim.Const(value =L_4)  # unit geometry 
sim.L_5 = sim.Const(value =L_5)  # unit geometry 
sim.dx  = sim.Const(value = dx)   # ft delta x

#Prime Mover Constants (Torque Balance)
sim.tau_p  = sim.Const(value = 3)   #tau
sim.k_gain = sim.Const(value = 1)      #one to one ratio between torque and SPM

##Economic
sim.Weight_lb_ft = sim.Const(value = sim.rho_r.value*sim.Ac.value*sim.g.value/sim.g_conv/144)    #Weight of rod string, lbf/ft
sim.Capex        = sim.Const(value = 200000)          #Cost of Pumping Rod Unit,$?
sim.P_o          = sim.Const(value = 50)               #Price of Oil, $/STB
sim.r            = sim.Const(value= .12/365)             #Daily Discount Rate, %
sim.P_th         = sim.Const(value = 100)             #tubing head pressure, psi
sim.TVD          = sim.Const(value = 4800)             #true vertical depth, ft              
sim.E_cost       = sim.Const(value = 0.13/3600)     #Cost of Electricity, cents/Kws

#Calculated Constants #DO NOT MODIFY#
sim.Wr          = sim.Const(value    = sim.TVD.value*sim.Weight_lb_ft.value)                                      #Weight of entire rod string, lbm
sim.D_a         = sim.Const(value   = 2*12*sim.rw.value)                                                    #Annulus Diameter, in
sim.gamma       = sim.Const(141.5/(sim.API.value+131.5))                                                #Specific gravity of Fluid
sim.P_startpump = sim.Const(value = 0.433*sim.gamma.value*sim.TVD.value)                              #Average Reservoir Pressure at Pump start up
sim.Pi          = sim.Const(value    = .433*sim.TVD.value)                                                    #Initial Reservoir Pressure, psi
sim.V_i         = sim.Const(value   = 7758*sim.A_d.value*sim.h_pz.value*sim.porosity.value*(1-sim.sw.value)/sim.Bo.value)     #OOIP, stb 


sim.A_t   = sim.Const((np.pi/4)*sim.D_t.value**2)  #Cross sectional Area of tubing, in^2
sim.A_a   = sim.Const((np.pi/4)*sim.D_a.value**2)  #Cross Sectional Area of Annulus, in^2
sim.Wf    = sim.Const(value = sim.TVD.value*sim.rho_w.value*sim.gamma.value*sim.g.value/sim.g_conv.value*(sim.A_t.value-sim.Ac.value)/144) # lbf, weight of fluid in tubing 

#MV's 
sim.SPM_in = sim.MV(value = 10, lb = 5, ub = 15)            #Rod Pump Pumping Speed/Torque, spm
                                               
#Variables
sim.Vp  = sim.Var(value = sim.V_i.value*(np.exp(sim.c.value*(sim.Pi.value-sim.P_startpump.value))-1))          #initial volume produced prior stb
sim.h   = sim.Var(value = 1.0*sim.TVD.value * 3 / 4)                                                        # Height, ft
sim.NPV = sim.Var(value = -1.0*sim.Capex.value)                                                   #Net Present Value, $
sim.y   = sim.Var( lb = -1, ub = 1) # SIGN(x)
sim.sa  = sim.Var(value = 0, lb = 0) # slack  variable a
sim.sb  = sim.Var(value = 0, lb = 0) # slack variable b
sim.tsi = sim.Var(value = 0.0) # simulation time
sim.SPM = sim.Var(value = 10)     #SPM, strokes/min
#omega = m.Var(value = 0)
sim.theta = sim.Var(value = 0) # rad i.e sec^-1 crank angle of surface unit
sim.u     = [sim.Var(value = 9.22,name  = 'un'+str(i)) for i in range(npx) ] # relative position of each rod segment
sim.v     = [sim.Var(value = 0.0, name  = 'vn'+str(i)) for i in range(npx)] # velocity of reach rod segment
sim.f     = [sim.Var(value = 0.0, name  = 'fn'+str(i)) for i in range (npx)] # load at each rod segment
sim.P     = sim.Var(value = 1e-6) # unitless, load at the pump 


## State Variables
sim.P_res = sim.Var(value = sim.P_startpump.value*1.0)                                                                                                                        #Current Reservoir Pressure , psi
sim.P_wf  = sim.Var(value = 0.433*sim.gamma*sim.h.value)                                                                                                                         #Bottomhole Flowing Pressure, psi
sim.q_in  = sim.Var(value = (1/86400)*sim.k.value*sim.h_pz.value*(sim.P_res.value-sim.P_wf.value)/(141.2*sim.Bo.value*sim.u_visc.value*((1/2)*np.log(4*sim.A_d.value/(sim.gamma_E.value*sim.C_a.value*sim.rw.value**2)) + sim.S.value))) #IPR-VLP Flow rate, STB/s
sim.q_out = sim.Var(value = 0)        # Outgoing Flow Rate, STB/s
sim.t     = sim.Var(value = 0)                                                                                                                                                #Time, days
sim.W_rod = sim.Var(value = (1.0962)*sim.q_out.value*(sim.P_th.value-sim.P_wf.value + .433*sim.gamma.value*sim.TVD.value) + (4.7053e-7)*sim.Wr.value*sim.St_length.value*sim.SPM.value)                    #Work supplied by electric Motor, KW

#Intermediates
sim.hs   = sim.Intermediate(sim.sqrt(L_1**2 +L_2**2 + 2 *L_1 *L_2 *sim.cos(sim.theta)))
sim.V_rp = sim.Intermediate((1/9702)*(np.pi/4)*sim.D_t**2*sim.St_length)   #Volume Extracted per stroke length, STB

#Equations
##AlgebraicEqns 
sim.Equation(sim.P_wf  == 0.433*sim.gamma*sim.h)
sim.Equation(sim.P_res == sim.Pi-(1/sim.c)*sim.log((sim.Vp/sim.V_i)+1))
sim.Equation(sim.q_in  == (1/86400)*sim.k*sim.h_pz*(sim.P_res-sim.P_wf)/(141.2*sim.Bo*sim.u_visc*((1/2)*sim.log(4*sim.A_d/(sim.gamma_E*sim.C_a*sim.rw**2)) + sim.S)))  #STB/s
sim.Equation(sim.W_rod == (1.0962)*sim.q_out*(sim.P_th-sim.P_wf + .433*sim.gamma*sim.TVD) + (4.7053e-7)*sim.Wr*sim.St_length*sim.SPM)  

#Prime Mover Equations- Torque Balance and Kinematic Eqns
sim.Equation(sim.SPM.dt()          == -(1/sim.tau_p)*sim.SPM + (sim.k_gain/sim.tau_p)*sim.SPM_in)
sim.Equation((2*sim.pi/60)*sim.SPM == sim.theta.dt())

#m.Equation(theta ==2.5118541087922712 +tsi*SPM_in * 2.0 * 3.147 / 60.0) # convert time to angle in radians
#m.Equation(SPM == omega/(2*pi)/60)

sim.Equation(sim.u[0]        == (1/12)*L_5*(sim.asin(L_1*sim.sin(sim.theta)/sim.hs)+sim.acos((sim.hs**2+L_3**2-L_4**2)/(2*L_3*sim.hs)))) # position of polished rod, inches 
[sim.Equation(sim.v[i+1].dt()== sim.a**2 * (sim.u[i+2] - 2.0*sim.u[i+1] + sim.u[i])/sim.dx**2 - sim.pi*sim.a*sim.nu/(2.0*sim.TVD)*sim.v[i+1] - (1-sim.rho_w*sim.gamma/sim.rho_r)*sim.g) for i in range(npx-2) ]# wave equation
sim.Equation(sim.q_out       == sim.A_t * sim.u[-1].dt()*12/231/42 * (1+sim.y)/2) # rate of fluid production, barrels/
#m.Equation(q_out == (1/60)*V_rp*SPM)

# Equations for calculating rod loading
# Load at surface
sim.Equation(sim.f[0] == sim.E*sim.Ac*1/2/sim.dx *(-sim.u[2] + 4*sim.u[1] -3*sim.u[0]))
# Load at pump
#m.Equation(f[npx-1] == E*Ac* 1/2.0/dx *(3*u[npx-1] - 4*u[npx-2] + u[npx-3]))
sim.Equation(sim.f[npx-1] == sim.E*sim.Ac* sim.P)
# load at intermediate points
[sim.Equation(sim.f[1+i] == sim.E*sim.Ac*1/2.0/dx*(sim.u[i+2] - sim.u[i])) for i in range(npx-2)]
# pump boundary
sim.Equation( sim.u[npx-1]*sim.alpha + (sim.u[npx-1] - sim.u[npx-2])/dx == sim.P)
#add in signum for lifting and lowering conditions
sim.Equation(sim.v[-1] == sim.sb - sim.sa )
sim.Equation(sim.P == -((sim.Wf- (sim.A_t - sim.Ac)*sim.P_wf)/sim.E/sim.Ac) * (1 + sim.y)/2 ) # -P_wf*A_t


##DifferentialEans
sim.Equation(sim.t.dt() == 1)
sim.Equation(sim.Vp.dt() == sim.q_in)
sim.Equation(sim.NPV.dt() == (sim.P_o*sim.q_out-sim.E_cost*sim.W_rod)*sim.exp(-sim.r*sim.t))
sim.Equation(sim.h.dt() == (1617/2)*(sim.q_in - sim.q_out)/(sim.A_a -sim.A_t))
sim.Equation(sim.tsi.dt()==1.0) # create time variable
[sim.Equation(sim.u[i].dt()==sim.v[i]) for i in range(npx)] # velocity of rod string 

# Set Objectives ##################################################
sim.Obj((sim.sa*(1+sim.y) + sim.sb*(1-sim.y))) # objective function to make signum work.
sim.Equation((sim.sa*(1+sim.y) + sim.sb*(1-sim.y))<= 1e-6)

#%%
#SetGlobalOptions Simulator ##############################################################
sim.options.IMODE = 5     # 4 = Dynamic Simulation (Seqential)
sim.options.NODES = 2    # 3 = 3 Nodes, 2 = No collocation nodes
sim.options.SOLVER = 3    # 1 =APOPT, 3 = IPOPT
sim.options.time_shift = npt-1 # time shift forward for multiple simulations
sim.options.MAX_ITER = 450

#SetLocalOptions Simulator###############################################################
sim.SPM_in.FSTATUS = 0#1 # accept measurments 
sim.SPM_in.STATUS  = 0 # don't let optimizer change (simulation)

#%%
# Simulate the application in loop

loops = 1800 # number of steps forward in time
sim_time = tf * loops # Total length of simulation horizon (sec)

time_loops1800 = np.zeros(loops) # store the time to do each solution

plt.figure()
plt.ion()
plt.show()

for i in range(loops):
    
    # start loop time
    start_loop_time = time.time()
    
    #############################

    #Doublet Test ###############
    
    if i > loops*(3/4):
        sim.SPM_in.VALUE = 10
    
    if i > loops*(2/4) and i < loops*(3/4):
        sim.SPM_in.VALUE = 5
        
    if i > (1/4)*loops and i < loops*(2/4):
        sim.SPM_in.VALUE = 15
    
    if i < (1/4)*loops or i < 10:#10: # TODO: if rerunnign this probably change this value of 61 back to 10
        sim.SPM_in.VALUE = 10#10
        
    #############################

    # simulate system for 1 second
    sim.solve()
    if i == 0:
        # Create and store results
        ts = np.array(sim.tsi.value) # simulation time storage
        us = [np.array(sim.u[i].value) for i in range(npx)] # u relative position storage
        vs = [np.array(sim.v[i].value) for i in range(npx)]
        fs = [np.array(sim.f[i].value) for i in range(npx)] # dynamic load storage 
        hstor = np.array(sim.h.value) # height of fluid in annulus storage
        q_ins= np.array(sim.q_in.value) # reservoir influx storage
        q_outs = np.array(sim.q_out.value) # production rate storage
        P_ress = np.array(sim.P_res.value) # reservoir pressure storage
        Vps    = np.array(sim.Vp.value) # cumulative volume produced storage
        NPVs   = np.array(sim.NPV.value) # NPV storage
        W_rods = np.array(sim.W_rod.value) # work of rod (work to lift fluid) storage
        SPMs   = np.array(sim.SPM_in.value) # Strokes per minute/ Torque storage Set Points
        SPMr   = np.array(sim.SPM.value)    #SPM storage
        thetas   = np.array(sim.theta.value)#Theta storage
        P_wfs  = np.array(sim.P_wf.value) # bottom hole pressure storage
        ys     = np.array(sim.y.value) # sign of du/dt storage

    elif i>0:
        ts = np.append(ts,sim.tsi.value) # simulation time storage
        us = [np.append(us[i],sim.u[i].value) for i in range(npx)] # u relative position storage
        vs = [np.append(vs[i],sim.v[i].value) for i in range(npx)]
        fs = [np.append(fs[i],sim.f[i].value) for i in range(npx)] # dynamic load storage 
        hstor = np.append(hstor,sim.h.value) # height of fluid in annulus storage
        q_ins= np.append(q_ins,sim.q_in.value) # reservoir influx storage
        q_outs = np.append(q_outs,sim.q_out.value) # production rate storage
        P_ress = np.append(P_ress,sim.P_res.value) # reservoir pressure storage
        Vps    = np.append(Vps,sim.Vp.value) # cumulative volume produced storage
        NPVs   = np.append(NPVs,sim.NPV.value) # NPV storage
        W_rods = np.append(W_rods,sim.W_rod.value) # work of rod (work to lift fluid) storage
        SPMs   = np.append(SPMs,sim.SPM_in.value) # Strokes per minute storage
        SPMr   = np.append(SPMr,sim.SPM.value)      #Strokes per minute storage
        thetas = np.append(thetas,sim.theta.value)  
        P_wfs  = np.append(P_wfs,sim.P_wf.value) # bottom hole pressure storage
        ys     = np.append(ys,sim.y.value) # sign of du/dt storage
    
    # Plotting
    plt.clf()
    ax=plt.subplot(311)
    ax.grid()
    plt.plot(ts[0:i*npt],SPMs[0:i*npt],'ro',label='Motor Torque')
    plt.plot(ts[0:i*npt],SPMr[0:i*npt],'bo',label='SPM')
    plt.ylabel('Strokes per Minute')
    plt.legend(loc=2)
    ax=plt.subplot(312)
    ax.grid()
    plt.plot(ts[0:i*npt],hstor[0:i*npt],'k-',label= 'height')
    plt.plot(ts[0:i*npt], np.ones(i*npt)*(sim.TVD.value*3/4 -3), label = 'height SP')
    plt.ylabel('Annular Fluid Height')
    plt.legend(loc='best')
    ax = plt.subplot(313)
    ax.grid()
    plt.plot(ts[0:i*npt], q_outs[0:i*npt], label = 'q_out')
    plt.plot(ts[0:i*npt], q_ins[0:i*npt], label = 'q_in')
    plt.legend()
    plt.ylabel('Flow Rate, STB/s')
    plt.xlabel('Time (sec)')
    plt.draw()
    plt.pause(0.02)
###############################################################
    # end loop time
    end_loop_time = time.time()
    time_loops1800[i] = end_loop_time - start_loop_time

#PlotResults###################################################################
#%%
## Figure 1
plt.figure(figsize = (6,4.5))
plt.subplot(211)
plt.plot(ts, hstor, 'r--', label = 'height in annulus')
plt.ylabel('height, ft')
plt.legend()

plt.subplot(212)
plt.plot(ts, q_ins, 'b--', label = r'$q_{in}$')
plt.plot(ts, q_outs, 'g--', label = r'$q_{out}$')
plt.ylabel('Flow Rate, STB/s')
plt.xlabel('time, sec')
plt.legend()

plt.show()

##Figure 2
plt.figure(figsize = (6,4.5))
plt.subplot(211)
plt.plot(ts, P_ress, 'k--', label = 'Reservoir Pressure')
#plt.plot(m.time, P_wf.value, 'r--', label = r'$P_{wf}$')
plt.ylabel('Pressure, psi')
plt.legend()

plt.subplot(212)
plt.plot(ts, Vps, '--', label = 'Cumulative Volume Produced')
plt.ylabel('Volume, STB')
plt.xlabel('time, sec')
plt.legend()
plt.tight_layout()

plt.show()

##Figure 3
plt.figure(figsize = (6,4.5))
plt.plot(ts, NPVs/(1e6), 'g:', label = 'NPV')
plt.xlabel('time, sec')
plt.ylabel('NPV, $ Millions')
plt.legend()

plt.show()

#Figure 4
plt.figure(figsize = (6,4.5))
plt.subplot(311)
plt.plot(ts,W_rods, 'b-', label = 'Work Supplied by Motor' )
plt.ylabel('KiloWatts, KW')

plt.subplot(312)
plt.plot(ts, SPMs, 'r-', label = 'Work Supplied by Motor' )
plt.ylabel('SPM')

plt.subplot(313)
plt.plot(ts, P_wfs, 'r--', label = r'$P_{wf}$')
plt.ylabel('FBHP, psi')
plt.xlabel('time, sec')
#plt.tight_layout()

plt.show()

##Figure 5 -Doublet Test
plt.figure(figsize = (6,4.5))
plt.subplot(211)
plt.plot(ts, hstor, 'r--', label = 'height in annulus')
plt.ylabel('height, ft')
plt.legend()

plt.subplot(212)
plt.plot(ts, SPMs, 'b--', label = r'SPM')
plt.ylabel('strokes/min')
plt.xlabel('time, sec')
plt.legend()

plt.show()

# store results in to structure for 3-d plotting
for i in range(npx):
    if i ==0:
        ustor = np.array([us[i]])
        tstor = np.array([ts])
    else:
        ustor = np.vstack([ustor,us[i]])
        tstor = np.vstack([tstor,ts])

for i in range(len(ts)):
    if i == 0:
        xstor = xi
    else:
        xstor = np.vstack([xstor,xi])
x = xstor.T
t = tstor
ustor = np.array(ustor)

fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection='3d')
ax.set_xlabel('Distance x')
ax.set_ylabel('Time t')
p = ax.plot_wireframe(x,t,ustor,rstride=1,cstride=1)
fig.show()  

plt.figure()
[plt.scatter(ts,us[i], label ='u' +str(i)) for i in range(npx)]
#plt.plot(m.time,u[-1],'r--')
plt.legend()
plt.show()

# Plot surface dynagraph ()
plt.figure()
plt.plot((us[0]- np.min(us[0]))*12,-fs[0],label = 'Surface Dynagraph') # dynamic plus static load #+ sim.TVD.value*sim.Ac.value*sim.rho_r.value/144
plt.legend()
plt.xlabel('Position (in)')
plt.ylabel('Load (lbf)')
plt.show()

# plot pump dynagraph
plt.figure()
plt.plot((us[npx-1]-np.min(us[npx-1]))*12,-fs[npx-1], label = 'Pump Dynagraph')
plt.xlabel('Position (in)')
plt.ylabel('Load (lbf)')
plt.legend()
plt.show()

# plot pump position vs. time 
plt.figure()
plt.plot(ts,us[-1], label = 'Pump Position')
plt.plot(ts,vs[-1], label ='Pump Velocity')
plt.plot(ts,np.array(ys)*10, label = 'Sign(V)')
plt.legend()
plt.show()

plt.figure()
plt.plot(ts,-fs[-1], label = 'Pump Load (lbf)')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Load (lbf)')
plt.legend()
plt.tight_layout()

plt.show()

#Plot SPM Dynamics
plt.figure()
plt.plot(ts,SPMs)
plt.plot(ts,SPMr)
plt.ylabel('SPM')
plt.xlabel('time,s')

plt.show()

#Plot theta Dynamics (radians)
plt.figure()
plt.plot(ts,thetas)
plt.xlabel('time, s')
plt.ylabel('radians')
plt.show()

# endtime
end_time = time.time()
total_time1800 = end_time - start_time

#%% # updated graphs for putting in the paper

# Plot surface dynagraph () where it works
plt.figure()
plt.plot((us[0][-600:]- np.min(us[0][-600:]))*12,-fs[0][-600:],label = 'Surface Dynagraph') # dynamic plus static load #+ sim.TVD.value*sim.Ac.value*sim.rho_r.value/144
plt.legend()
plt.xlabel('Position (in)')
plt.ylabel('Load (lbf)')
plt.show()

# plot pump dynagraph where it works
plt.figure()
plt.plot((us[npx-1][-600:]-np.min(us[npx-1][-600:]))*12,-fs[npx-1][-600:], label = 'Pump Dynagraph')
plt.xlabel('Position (in)')
plt.ylabel('Load (lbf)')
plt.legend()
plt.show()

# plot pump position vs. time 
plt.figure()
plt.plot(ts,us[-1], label = 'Pump Position')
plt.plot(ts,vs[-1], label ='Pump Velocity')
plt.plot(ts,np.array(ys)*10, label = 'Sign(V)')
plt.legend()
plt.show()

# time for pump dynagraph
plt.figure()
plt.plot(ts,-fs[-1], label = 'Pump Load (lbf)')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Load (lbf)')
plt.legend()
plt.tight_layout()
plt.show()

# time for surface dynagraph
plt.figure()
plt.plot(ts,-fs[0], label = 'Pump Load (lbf)')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Load (lbf)')
plt.legend()
plt.tight_layout()
plt.show()

#%% Store Results in a dictionary
    
res['ts'        ] = ts
res['us'        ] = us
res['vs'        ] = vs
res['fs'        ] = fs
res['hstor'     ] = hstor
res['q_ins'     ] = q_ins
res['q_outs'    ] = q_outs
res['P_ress'    ] = P_ress
res['Vps'       ] = Vps
res['NPVs'      ] = NPVs
res['W_rods'    ] = W_rods   
res['SPMs'      ] = SPMs
res['SPMr'      ] = SPMr
res['thetas'    ] = thetas
res['P_wfs'     ] = P_wfs
res['ys'        ] = ys
res['times'     ] = time_loops1800
res['total_time'] = total_time1800

np.save('Simulation_data_' + str(loops) + 's.npy', res)

#%% Load results from the dictionary

res = np.load('Simulation_data_1800s.npy').item() # res = np.load('Simulation_data_1800s.npy').item()
loops = 1800

#Define Timespace
tf  = 1 # sec, length of simulation 
npt = 30 # number of time discretizations
nit = 1 # number of itterations to solve 
ti = np.linspace(0,tf,npt) # times for plotting
                        
# Define Rod Discretizations
TVD_d = 4800 # ft, lenth of rod
npx = 30 # number of rod discretizations
dx = TVD_d/(npx-1) # #ft lenth of rod discretizations
xi = np.linspace(0,TVD_d,npx) # possitions allong rod (for plotting)

fig = plt.figure() # see line 673

#%% plot using the dictionary
# figure 1

plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.formatter.useoffset'] = False
plt.subplot(211)
plt.plot(res['ts'], res['hstor'], 'r--', label = 'Fluid Level (ft)')
plt.ylabel('Height (ft)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)

plt.subplot(212)
plt.plot(res['ts'], res['q_ins'], 'b--', label = r'$q_{in}$')
plt.plot(res['ts'], res['q_outs'], 'g--', label = r'$q_{out}$')
plt.ylabel('Flow (STB/s)', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.tight_layout()

plt.savefig('fig_1.eps', dpi = 1200, Transparent = True)
plt.show()

#Figure 2

plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.subplot(211)
plt.plot(res['ts'], res['P_ress'], 'k--', label = 'Reservoir Pressure')
#plt.plot(m.time, P_wf.value, 'r--', label = r'$P_{wf}$')
plt.ylabel('Pressure (psi)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)

plt.subplot(212)
plt.plot(res['ts'], res['Vps'], '--', label = 'Cumulative Volume Produced')
plt.ylabel('Volume (STB)', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.tight_layout()

plt.savefig('fig_2.eps', dpi = 1200, Transparent = True)
plt.show()

#Figure 3

plt.figure()
plt.plot(res['ts'], res['NPVs']/(1e6), 'g:', label = 'NPV')
plt.xlabel('Time (seconds)', fontsize = 12)
plt.ylabel('NPV ($ Millions)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.tight_layout()
plt.savefig('fig_3.eps', dpi = 1200, Transparent = True)
plt.show()

#%% Figure 4
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.formatter.useoffset'] = False
plt.subplot(311)
plt.plot(res['ts'],res['W_rods'], 'b-', label = 'Lifting Power' )
plt.ylabel('KiloWatts (kW)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.xlim(0,loops)

plt.subplot(312)
plt.plot(res['ts'], res['SPMs'], 'r--', label = r'$T_{net}$')
# add Tnet
plt.plot(res['ts'], res['SPMr'], 'b-', label = r'$SPM$')
#
plt.ylabel('SPM', fontsize = 12)
plt.legend(loc='upper right', fontsize = 12)
plt.xlim(0,loops)

plt.subplot(313)
plt.plot(res['ts'], res['P_wfs'], 'r--', label = r'$P_{wf}$')
plt.ylabel('FBHP (psi)', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.xlim(0,loops)
plt.tight_layout()

plt.savefig('fig_4.eps', dpi = 1200, Transparent = True)
plt.show()
#%%

# Figure 5 -Doublet Test

plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.formatter.useoffset'] = False
plt.subplot(211)
plt.plot(res['ts'], res['hstor'], 'r--', label = 'Height in Annulus')
plt.ylabel('Height (ft)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)

plt.subplot(212)
plt.plot(res['ts'], res['SPMs'], 'b--', label = r'$T_{net}$')
plt.ylabel(r'$T_{net} (ft-lbs)$', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)

plt.tight_layout()
plt.savefig('fig_5.eps', dpi = 1200, Transparent = True)
plt.show()

# Figure 6

# store results in to structure for 3-d plotting
for i in range(npx):
    if i ==0:
        ustor = np.array([res['us'][i]]) # TODO: check and fix?
        tstor = np.array([res['ts']])
    else:
        ustor = np.vstack([ustor,res['us'][i]])
        tstor = np.vstack([tstor,res['ts']])

for i in range(len(res['ts'])):
    if i == 0:
        xstor = xi
    else:
        xstor = np.vstack([xstor,xi])
x = xstor.T
t = tstor
ustor = np.array(ustor)

fig # note it takes a long time to make this graph
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
ax = fig.add_subplot(1,1,1,projection='3d')
ax.set_xlabel('Distance (ft)')
ax.set_ylabel('Time (seconds)')
ax.set_zlabel('Position (ft)')
p = ax.plot_wireframe(x,t,ustor - np.min(ustor),rstride=1,cstride=1)
fig.savefig('fig_6.eps', dpi = 1200, Transparent = True)
fig.show()  

#figure 7

plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
[plt.scatter(res['ts'],res['us'][i], label ='u' + str(i)) for i in range(npx)]
#plt.plot(m.time,u[-1],'r--')
plt.legend(loc='best', fontsize = 12)
plt.savefig('fig_7.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 8

# Plot surface dynagraph ()
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot((res['us'][0]- np.min(res['us'][0]))*12,-res['fs'][0],label = 'Surface Dynagraph') # dynamic plus static load #+ sim.TVD.value*sim.Ac.value*sim.rho_r.value/144
plt.legend(loc='best', fontsize = 12)
plt.xlabel('Position (in)', fontsize = 12)
plt.ylabel('Load (lbf)', fontsize = 12)
plt.savefig('fig_8.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 9

# plot pump dynagraph
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot((res['us'][npx-1]-np.min(res['us'][npx-1]))*12,-res['fs'][npx-1], label = 'Pump Dynagraph')
plt.xlabel('Position (in)', fontsize = 12)
plt.ylabel('Load (lbf)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.savefig('fig_9.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 10

# plot pump position vs. time 
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot(res['ts'],res['us'][-1], label = 'Pump Position')
plt.plot(res['ts'],res['vs'][-1], label ='Pump Velocity')
plt.plot(res['ts'],np.array(res['ys'])*10, label = 'Sign(V)')
plt.legend(loc='best', fontsize = 12)
fig.savefig('fig_10.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 11

plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot(res['ts'],-res['fs'][-1], label = 'Pump Load (lbf)')
plt.legend(loc='best', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.ylabel('Load (lbf)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.tight_layout()

plt.savefig('fig_11.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 12

#Plot SPM Dynamics
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot(res['ts'],res['SPMs'])
plt.plot(res['ts'],res['SPMr'])
plt.ylabel('SPM', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)

plt.savefig('fig_12.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 13

#Plot theta Dynamics (Radians)
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot(res['ts'],res['thetas'])
plt.xlabel('Time (seconds)', fontsize = 12)
plt.ylabel('Radians', fontsize = 12)
plt.savefig('fig_13.eps', dpi = 1200, Transparent = True)
plt.show()

#%% updated graphs for putting in the paper

# figure 14

# Plot surface dynagraph () where it works
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot((res['us'][0][-600:]- np.min(res['us'][0][-600:]))*12,-res['fs'][0][-600:],label = 'Surface Dynagraph') # dynamic plus static load #+ sim.TVD.value*sim.Ac.value*sim.rho_r.value/144
plt.legend(loc='best', fontsize = 12)
plt.xlabel('Position (in)', fontsize = 12)
plt.ylabel('Load (lbf)', fontsize = 12)
plt.tight_layout()
plt.savefig('fig_14.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 15

# plot pump dynagraph where it works
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot((res['us'][npx-1][-600:]-np.min(res['us'][npx-1][-600:]))*12,-res['fs'][npx-1][-600:], label = 'Pump Dynagraph')
plt.xlabel('Position (in)', fontsize = 12)
plt.ylabel('Load (lbf)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.tight_layout()
plt.savefig('fig_15.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 16

# plot pump position vs. time 
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot(res['ts'],res['us'][-1]-np.min(res['us'][-1]), label = 'Pump Position (ft)')
# the minus np.min(res['us'][-1] is used to put the minimum at zero
plt.plot(res['ts'],res['vs'][-1], label ='Pump Velocity (ft/s)')
plt.plot(res['ts'],np.array(res['ys'])*10, label = 'Sgn(V)')
plt.legend(loc='best', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.tight_layout()
plt.savefig('fig_16.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 17

# time for pump dynagraph
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot(res['ts'],-res['fs'][-1], label = 'Pump Load (lbf)')
plt.legend(loc='best', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.ylabel('Load (lbf)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.tight_layout()
plt.savefig('fig_17.eps', dpi = 1200, Transparent = True)
plt.show()

# figure 18

# time for surface dynagraph
plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.plot(res['ts'],-res['fs'][0], label = 'Pump Load (lbf)')
plt.legend(loc='best', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.ylabel('Load (lbf)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.tight_layout()
plt.savefig('fig_18.eps', dpi = 1200, Transparent = True)
plt.show()
