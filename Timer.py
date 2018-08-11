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
#Define Timespace for simulation
tf  = 1 # sec, length of simulation 
npt = 21 # number of time discretizations
nit = 1 # number of itterations to solve 
ti = np.linspace(0,tf,npt) # times for plotting

res = {} #set up on empty dictionary

# Define Rod Discretizations
TVD_d = 4800 # ft, lenth of rod
npx = 10 # number of rod discretizations
dx = TVD_d/(npx-1) # #ft lenth of rod discretizations
xi = np.linspace(0,TVD_d,npx) # possitions allong rod (for plotting)

#Set Points
SPM_in = np.ones(npt)*10

#BuildModle####################################################################

sim = GEKKO()

#Horizon Window  
sim.time = np.linspace(0,tf,npt)

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

starting_height = 3/4
#Simulation########################################

for m in [sim]:   
   #Constants
    m.API       = m.Const(value = 45)               #API gravity of fluid, unitless
    m.c         = m.Const(value = 0.000013)           #Compressibility, psi^-1
    m.k         = m.Const(value = 15)                 #Permeability, md
    m.Bo        = m.Const(value = 1.2)               #FVF, rb/STB
    m.A_d       = m.FV(value = 2, ub = 8, lb = 1)                #Drainage Area, Acres
    m.sw        = m.Const(value = 0.2)               #Water Saturation
    m.porosity  = m.FV(value = 0.08, ub = .12, lb = 0.07)        #Porosity, unitless
    m.gamma_E   = m.Const(value = 1.78)         #Euler Constant
    m.C_a       = m.Const(value = 31.6)             #Drainage Area Shape Factor (Circular)
    m.rw        = m.Const(value = 0.328)             #Welbore radius, ft
    m.S         = m.FV(value = 0, ub = 10, lb = -5)                  #unitless
    m.u_visc    = m.Const(value = 1.5)                # Viscosity, cp
    m.h_pz      = m.Const(value = 8)              #pay zone thickness, ft
    m.D_t       = m.Const(value = 2.5)              # tubing diameter, in
    m.St_length = m.Const(value = 85)         # rod pump stroke length, in
    
    m.g     = m.Const(value = 32.2)             # acceleration due to gravity, ft/s^3
    m.g_conv= m.Const(value = 32.2)           # lbf conversion , lb-ft/s^2-lbf
    m.rho_r = m.Const(value = 490)              # lbs/ft^3, density of rod steel
    m.rho_w = m.Const(value = 62.3 )          # lbs/ft^3, density of water at standard conditions
    m.a     = m.Const(value =18996.06 )          # ft/s speed of sound in steel
    m.D_r   = m.Const(value = 1.0)               # in, diameter of rod string
    m.Ac    = m.Const(value= m.D_r.value**2/4.0*np.pi) # in^2, cross sectional area of rod 
    m.nu    = m.Const(value = 0.01)              # unitless, damping coefficient
    m.pi    = m.Const(value=np.pi)
    m.E     = m.Const(value = 32025000.0)        # psi sucker rod modulus of elasticity
    m.alpha = m.Const(value = 0.0)            # pump parameter, unitless
    m.beta  = m.Const(value = 1.0)            # pump parameter, unitless 
    
    m.L_1 = m.Const(value =L_1) # unit geometry 
    m.L_2 = m.Const(value =L_2)  # unit geometry 
    m.L_3 = m.Const(value =L_3)  # unit geometry 
    m.L_4 = m.Const(value =L_4)  # unit geometry 
    m.L_5 = m.Const(value =L_5)  # unit geometry 
    m.dx  = m.Const(value = dx)   # ft delta x
    
    #Prime Mover Constants (Torque Balance)
    m.tau_p = m.Const(value = 3)   #tau
    m.k_gain = m.Const(value = 1)      #one to one ratio between torque and SPM
    
    ##Economic
    m.Weight_lb_ft = m.Const(value = m.rho_r.value*m.Ac.value*m.g.value/m.g_conv/144)    #Weight of rod string, lbf/ft
    m.Capex        = m.Const(value = 200000)          #Cost of Pumping Rod Unit,$?
    m.P_o          = m.Const(value = 50)               #Price of Oil, $/STB
    m.r            = m.Const(value= .12/365)             #Daily Discount Rate, %
    m.P_th         = m.Const(value = 100)             #tubing head pressure, psi
    m.TVD          = m.Const(value = 4800)             #true vertical depth, ft              
    m.E_cost       = m.Const(value = 0.13/3600)     #Cost of Electricity, cents/Kws
    
    #Calculated Constants #DO NOT MODIFY#
    m.Wr          = m.Const(value    = m.TVD.value*m.Weight_lb_ft.value)                                      #Weight of entire rod string, lbm
    m.D_a         = m.Const(value   = 2*12*m.rw.value)                                                    #Annulus Diameter, in
    m.gamma       = m.Const(141.5/(m.API.value+131.5))                                                #Specific gravity of Fluid
    m.P_startpump = m.Const(value = 0.433*m.gamma.value*m.TVD.value)                              #Average Reservoir Pressure at Pump start up
    m.Pi          = m.Const(value    = .433*m.TVD.value)                                                    #Initial Reservoir Pressure, psi
    
    
    m.A_t = m.Const((np.pi/4)*m.D_t.value**2)  #Cross sectional Area of tubing, in^2
    m.A_a = m.Const((np.pi/4)*m.D_a.value**2)  #Cross Sectional Area of Annulus, in^2
    m.Wf    = m.Const(value = m.TVD.value*m.rho_w.value*m.gamma.value*m.g.value/m.g_conv.value*(m.A_t.value-m.Ac.value)/144) # lbf, weight of fluid in tubing

    #MV's 
    m.SPM_in = m.MV(value = 15, lb = 0, ub = 15)            #Rod Pump Pumping Speed/Torque, spm
                                                   
    #Variables
    m.V_i= m.Var(value   = 7758*m.A_d.value*m.h_pz.value*m.porosity.value*(1-m.sw.value)/m.Bo.value)    #OOIP, stb
    m.Vp = m.Var(value = m.V_i.value*(np.exp(m.c.value*(m.Pi.value-m.P_startpump.value))-1))          #initial volume produced prior stb
    
    
    m.h = m.CV(value = 1.0*m.TVD.value*starting_height) # Height, ft
                                                             
    
    m.NPV = m.Var(value = -1.0*m.Capex.value)                                                   #Net Present Value, $
    m.y = m.Var( lb = -1, ub = 1) # SIGN(x)
    m.sa = m.Var(value = 0, lb = 0) # slack  variable a
    m.sb = m.Var(value = 0, lb = 0) # slack variable b
    m.tsi = m.Var(value = 0.0) # mulation time
    m.SPM = m.Var(value = 15)     #SPM, strokes/min
    #omega = m.Var(value = 0)
    m.theta = m.Var(value = 0) # rad i.e sec^-1 crank angle of surface unit
    m.u = [m.SV(value = 9.22) for i in range(npx)] # relative position of each rod segment
    m.v = [m.Var(value = 0.0) for i in range(npx)] # velocity of reach rod segment
    m.f = [m.SV(value = 0.0) for i in range (npx)] # load at each rod segment
    m.P = m.Var(value = 1e-6) # unitless, load at the pump 
    
    
    ## State Variables
    m.P_res = m.Var(value = m.P_startpump.value*1.0)                                                                                                                        #Current Reservoir Pressure , psi
    m.P_wf = m.Var(value = 0.433*m.gamma*m.h.value)                                                                                                                         #Bottomhole Flowing Pressure, psi
    m.q_in = m.Var(value = (1/86400)*m.k.value*m.h_pz.value*(m.P_res.value-m.P_wf.value)/(141.2*m.Bo.value*m.u_visc.value*((1/2)*np.log(4*m.A_d.value/(m.gamma_E.value*m.C_a.value*m.rw.value**2)) + m.S.value))) #IPR-VLP Flow rate, STB/s
    m.q_out = m.Var(value = 0)        # Outgoing Flow Rate, STB/s
    m.t = m.Var(value = 0)                                                                                                                                                #Time, days
    m.W_rod = m.Var(value = (1.0962)*m.q_out.value*(m.P_th.value-m.P_wf.value + .433*m.gamma.value*m.TVD.value) + (4.7053e-7)*m.Wr.value*m.St_length.value*m.SPM.value)                    #Work supplied by electric Motor, KW
    
    #Intermediates
    m.hs = m.Intermediate(m.sqrt(L_1**2 +L_2**2 + 2 *L_1 *L_2 *m.cos(m.theta)))
    m.V_rp = m.Intermediate((1/9702)*(np.pi/4)*m.D_t**2*m.St_length)   #Volume Extracted per stroke length, STB
    
    #Equations
    ##AlgebraicEqns
    m.Equation(m.V_i == 7758*m.A_d*m.h_pz*m.porosity*(1-m.sw)/m.Bo)
    m.Equation(m.P_wf == 0.433*m.gamma*m.h)
    m.Equation(m.P_res == m.Pi-(1/m.c)*m.log((m.Vp/m.V_i)+1))
    m.Equation(m.q_in == (1/86400)*m.k*m.h_pz*(m.P_res-m.P_wf)/(141.2*m.Bo*m.u_visc*((1/2)*m.log(4*m.A_d/(m.gamma_E*m.C_a*m.rw**2)) + m.S)))  #STB/s
    m.Equation(m.W_rod == (1.0962)*m.q_out*(m.P_th-m.P_wf + .433*m.gamma*m.TVD) + (4.7053e-7)*m.Wr*m.St_length*m.SPM)  
    
    #Prime Mover Equations- Torque Balance and Kinematic Eqns
    m.Equation(m.SPM.dt() == -(1/m.tau_p)*m.SPM + (m.k_gain/m.tau_p)*m.SPM_in)
    m.Equation((2*m.pi/60)*m.SPM == m.theta.dt())
    
    #m.Equation(theta ==2.5118541087922712 +tsi*SPM_in * 2.0 * 3.147 / 60.0) # convert time to angle in radians
    #m.Equation(SPM == omega/(2*pi)/60)
    
    m.Equation(m.u[0]        == (1/12)*L_5*(m.asin(L_1*m.sin(m.theta)/m.hs)+m.acos((m.hs**2+L_3**2-L_4**2)/(2*L_3*m.hs)))) # position of polished rod, inches 
    [m.Equation(m.v[i+1].dt()== m.a**2 * (m.u[i+2] - 2.0*m.u[i+1] + m.u[i])/m.dx**2 - m.pi*m.a*m.nu/(2.0*m.TVD)*m.v[i+1] - (1-m.rho_w*m.gamma/m.rho_r)*m.g) for i in range(npx-2) ]# wave equation
    m.Equation(m.q_out       == m.A_t * m.u[-1].dt()*12/231/42 * (1+m.y)/2) # rate of fluid production, barrels/
    #m.Equation(q_out == (1/60)*V_rp*SPM)
    
    # Equations for calculating rod loading
    # Load at surface
    m.Equation(m.f[0] == m.E*m.Ac*1/2/m.dx *(-m.u[2] + 4*m.u[1] -3*m.u[0]))
    # Load at pump
    #m.Equation(f[npx-1] == E*Ac* 1/2.0/dx *(3*u[npx-1] - 4*u[npx-2] + u[npx-3]))
    m.Equation(m.f[npx-1] == m.E*m.Ac* m.P)
    # load at intermediate points
    [m.Equation(m.f[1+i] == m.E*m.Ac*1/2.0/dx*(m.u[i+2] - m.u[i])) for i in range(npx-2)]
    # pump boundary
    m.Equation( m.u[npx-1]*m.alpha + (m.u[npx-1] - m.u[npx-2])/dx == m.P)
    #add in signum for lifting and lowering conditions
    m.Equation(m.v[-1] == m.sb - m.sa )
    m.Equation(m.P == -((m.Wf- (m.A_t - m.Ac)*m.P_wf)/m.E/m.Ac) * (1 + m.y)/2 ) # -P_wf*A_t
    
    
    ##DifferentialEans
    m.Equation(m.t.dt() == 1)
    m.Equation(m.Vp.dt() == m.q_in)
    m.Equation(m.NPV.dt() == (m.P_o*m.q_out-m.E_cost*m.W_rod)*m.exp(-m.r*m.t))
    m.Equation(m.h.dt() == (1617/2)*(m.q_in - m.q_out)/(m.A_a -m.A_t))
    m.Equation(m.tsi.dt()==1.0) # create time variable
    [m.Equation(m.u[i].dt()==m.v[i]) for i in range(npx)] # velocity of rod string 
    
    # Set Objectives ##################################################
    m.Obj((m.sa*(1+m.y) + m.sb*(1-m.y))) # objective function to make signum work.
    #m.Equation((m.sa*(1+m.y) + m.sb*(1-m.y))<= 1e-6)


#SetGlobalOptions(Simulation)##############################################################
sim.options.IMODE = 5     # 4 = Dynamic Simulation (Seqential)
sim.options.NODES = 2    # 3 = 3 Nodes, 2 = No collocation nodes
sim.options.SOLVER = 3    # 1 =APOPT, 3 = IPOPT
sim.options.time_shift = npt-1 # time shift forward for multiple simulations
sim.options.MAX_ITER = 450

#SetLocalOptions###############################################################
#N/A
sim.SPM_in.FSTATUS = 1 # accept measurments 
sim.SPM_in.STATUS  = 0 # don't let optimizer change (simulation)



#Solve#########################################################################
#%%
# Solve the simulation in a loop to simulate a longer horizon
loops = 180 # number of steps forward in time (previously 145 sec)




###############################################################
trigger = 30 # ammount of time to shut unit down when low level is detected
for i in range(loops):
    
   # simulate system for 1 second
    sim.solve()
    if i == 0:
        # Create and store results
        sim_ts = np.array(sim.tsi.value) # simulation time storage
        sim_us = [np.array(sim.u[i].value) for i in range(npx)] # u relative position storage
        sim_vs = [np.array(sim.v[i].value) for i in range(npx)]
        sim_fs = [np.array(sim.f[i].value) for i in range(npx)] # dynamic load storage 
        sim_hstor = np.array(sim.h.value) # height of fluid in annulus storage
        sim_q_ins= np.array(sim.q_in.value) # reservoir influx storage
        sim_q_outs = np.array(sim.q_out.value) # production rate storage
        sim_P_ress = np.array(sim.P_res.value) # reservoir pressure storage
        sim_Vps    = np.array(sim.Vp.value) # cumulative volume produced storage
        sim_NPVs   = np.array(sim.NPV.value) # NPV storage
        sim_W_rods = np.array(sim.W_rod.value) # work of rod (work to lift fluid) storage
        sim_SPMs   = np.array(sim.SPM_in.value) # Strokes per minute/ Torque storage Set Points
        sim_SPMr   = np.array(sim.SPM.value)    #SPM storage
        sim_thetas   = np.array(sim.theta.value)#Theta storage
        sim_P_wfs  = np.array(sim.P_wf.value) # bottom hole pressure storage
        sim_ys     = np.array(sim.y.value) # sign of du/dt storage

    elif i>0:
        sim_ts = np.append(sim_ts,sim.tsi.value) # simulation time storage
        sim_us = [np.append(sim_us[i],sim.u[i].value) for i in range(npx)] # u relative position storage
        sim_vs = [np.append(sim_vs[i],sim.v[i].value) for i in range(npx)]
        sim_fs = [np.append(sim_fs[i],sim.f[i].value) for i in range(npx)] # dynamic load storage 
        sim_hstor = np.append(sim_hstor,sim.h.value) # height of fluid in annulus storage
        sim_q_ins= np.append(sim_q_ins,sim.q_in.value) # reservoir influx storage
        sim_q_outs = np.append(sim_q_outs,sim.q_out.value) # production rate storage
        sim_P_ress = np.append(sim_P_ress,sim.P_res.value) # reservoir pressure storage
        sim_Vps    = np.append(sim_Vps,sim.Vp.value) # cumulative volume produced storage
        sim_NPVs   = np.append(sim_NPVs,sim.NPV.value) # NPV storage
        sim_W_rods = np.append(sim_W_rods,sim.W_rod.value) # work of rod (work to lift fluid) storage
        sim_SPMs   = np.append(sim_SPMs,sim.SPM_in.value) # Strokes per minute storage
        sim_SPMr   = np.append(sim_SPMr,sim.SPM.value)      #Strokes per minute storage
        sim_thetas = np.append(sim_thetas,sim.theta.value)  
        sim_P_wfs  = np.append(sim_P_wfs,sim.P_wf.value) # bottom hole pressure storage
        sim_ys     = np.append(sim_ys,sim.y.value) # sign of du/dt storage
    
    
    # ON/OFF control with timer
    
    h_meas = np.array(sim.h.value)
    h_min  = np.min(h_meas)
    
    trigger = trigger + 1
    if trigger >30:
        sim.SPM_in.value = 15
    
    if h_min < 3597:
        sim.SPM_in.value = 0.0
        trigger = 0
   
    if trigger<30:
             sim.SPM_in.value = 0.0
             
    #######################################################################
    
    # Plotting
    plt.clf()
    ax=plt.subplot(311)
    ax.grid()
    plt.plot(sim_ts[0:i*npt],sim_SPMs[0:i*npt],'ro',label='Motor Torque')
    plt.plot(sim_ts[0:i*npt],sim_SPMr[0:i*npt],'bo',label='SPM')
    plt.ylabel('Strokes per Minute')
    plt.legend(loc=2)
    ax=plt.subplot(312)
    ax.grid()
    plt.plot(sim_ts[0:i*npt],sim_hstor[0:i*npt],'k-',label= 'height')
    plt.plot(sim_ts[0:i*npt], np.ones(i*npt)*(sim.TVD.value*3/4 -3), label = 'height SP')
    #plt.plot(ts[0:i*npt], mpc_hs[0:npt*i], label = 'mpc height')
    plt.ylabel('Annular Fluid Height')
    plt.legend(loc='best')
    ax = plt.subplot(313)
    ax.grid()
    plt.plot(sim_ts[0:i*npt], sim_q_outs[0:i*npt], label = 'q_out')
    plt.plot(sim_ts[0:i*npt], sim_q_ins[0:i*npt], label = 'q_in')
    plt.legend()
    plt.ylabel('Flow Rate, STB/s')
    plt.xlabel('Time (sec)')
    plt.draw()
   
    plt.legend()
    plt.ylabel('Annular Fluid Height')
    plt.xlabel('Time (sec)')
    plt.draw()
    plt.pause(0.02)
    
# store in dictionary
#res['solve_stat' ] = solve_stat
res['sim_ts'      ] = sim_ts
res['sim_us'      ] = sim_us
res['sim_vs'      ] = sim_vs
res['sim_fs'      ] = sim_fs
res['sim_hstor'   ] = sim_hstor
res['sim_q_ins'  ] = sim_q_ins
res['sim_q_outs' ] = sim_q_outs
res['sim_P_ress'  ] = sim_P_ress
res['sim_Vps'     ] = sim_Vps
res[ 'NPVs'   ] = sim_NPVs
res['sim_W_rods'  ] = sim_W_rods   
res['sim_SPMs'    ] = sim_SPMs
res['sim_SPMr'    ] = sim_SPMr
res['sim_thetas'  ] = sim_thetas
res['sim_P_wfs'   ] = sim_P_wfs
res['sim_ys'      ] = sim_ys 

#%%

np.save('Timer_controller_' + str(loops) + 's21npt.npy', res)

#%%

sim_TVD_value = 4800 # this is a little weird and affects line 361 as well
res = np.load('Timer_controller_180s21npt.npy').item()

fig = plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

ax=plt.subplot(311)
ax.grid()
plt.plot(res['sim_ts'], res['sim_SPMs'], '--r', label=r'$T_{net}$ (ft-lb)') #r'Motor Torque ' for latex
plt.plot(res['sim_ts'], res['sim_SPMr'], '-b', label='SPM')
plt.ylabel('SPM', fontsize = 12)
plt.legend(loc=2, fontsize = 12)
plt.xlim(0,180)
ax=plt.subplot(312)
ax.grid()
plt.plot(res['sim_ts'], res['sim_hstor'],'k-',label= 'Actual')
plt.plot(res['sim_ts'], np.ones(np.size(res['sim_ts']))*(sim_TVD_value*3/4 -3), label = 'SP') # fix
#plt.plot(ts[0:i*npt], mpc_hs[0:npt*i], label = 'mpc height')
plt.ylabel('Fluid Level (ft)', fontsize = 12)
plt.legend(loc='best', fontsize = 12)
plt.xlim(0,180)
ax = plt.subplot(313)
ax.grid()
plt.plot(res['sim_ts'], res['sim_q_outs'], label = r'$q_{out}$')
plt.plot(res['sim_ts'], res['sim_q_ins'], label = r'$q_{in}$')
plt.legend(loc='best', fontsize = 12)
plt.xlim(0,180)
plt.ylabel('Flow (STB/s)', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.draw()
   
plt.legend(loc='best', fontsize = 12)
plt.xlim(0,180)
plt.ylabel('Flow (STB/s)', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)

plt.tight_layout()
fig.savefig('Timer_control.eps', dpi = 1200, Transparent = True)
plt.show()

##%% old graph
#
#plt.figure()
#
#ax=plt.subplot(311)
#ax.grid()
#plt.plot(sim_ts[0:i*npt],sim_SPMs[0:i*npt],'ro',label='Motor Torque')
#plt.plot(sim_ts[0:i*npt],sim_SPMr[0:i*npt],'bo',label='SPM')
#plt.ylabel('Strokes per Minute')
#plt.legend(loc=2)
#ax=plt.subplot(312)
#ax.grid()
#plt.plot(sim_ts[0:i*npt],sim_hstor[0:i*npt],'k-',label= 'height')
#plt.plot(sim_ts[0:i*npt], np.ones(i*npt)*(sim.TVD.value*3/4 -3), label = 'height SP')
##plt.plot(ts[0:i*npt], mpc_hs[0:npt*i], label = 'mpc height')
#plt.ylabel('Annular Fluid Height')
#plt.legend(loc='best')
#ax = plt.subplot(313)
#ax.grid()
#plt.plot(sim_ts[0:i*npt], sim_q_outs[0:i*npt], label = 'q_out')
#plt.plot(sim_ts[0:i*npt], sim_q_ins[0:i*npt], label = 'q_in')
#plt.legend()
#plt.ylabel('Flow Rate, STB/s')
#plt.xlabel('Time (sec)')
#plt.draw()
#   
#plt.legend()
#plt.ylabel('Annular Fluid Height')
#plt.xlabel('Time (sec)')
#
#plt.show()




