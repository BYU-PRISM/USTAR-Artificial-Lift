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
npt = 20 # number of time discretizations
nit = 1 # number of itterations to solve 
ti = np.linspace(0,tf,npt) # times for plotting


# Define time space for MPC
tf_mpc = 5
npt_mpc = npt*tf_mpc



# Define Rod Discretizations
TVD_d = 4800 # ft, lenth of rod
npx = 10 # number of rod discretizations
dx = TVD_d/(npx-1) # #ft lenth of rod discretizations
xi = np.linspace(0,TVD_d,npx) # possitions allong rod (for plotting)

#Set Points
SPM_in = np.ones(npt)*10


#BuildModle####################################################################
sim = GEKKO()
mhe = GEKKO()

#Horizon Window                             
sim.time = np.linspace(0,tf,npt)
mhe.time = np.linspace(0,tf,npt)

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

#Setpoints
level_height = 3
SP = starting_height*TVD_d-level_height
dSP = 0.3

#Simulation########################################
for m in [sim,mhe]:   
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
    m.SPM_in = m.MV(value = 15, lb = 5, ub = 15)            #Rod Pump Pumping Speed/Torque, spm
                                                   
    #Variables
    m.V_i= m.Var(value   = 7758*m.A_d.value*m.h_pz.value*m.porosity.value*(1-m.sw.value)/m.Bo.value)    #OOIP, stb
    m.Vp = m.Var(value = m.V_i.value*(np.exp(m.c.value*(m.Pi.value-m.P_startpump.value))-1))          #initial volume produced prior stb
    
    if m == sim:
        m.h = m.CV(value = 1.0*m.TVD.value*starting_height)
    else:
        m.h = m.MV(value = 1.0*m.TVD.value*starting_height, lb = 0, ub = 4800)                                                         # Height, ft
    
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
    
    m.Equation(m.u[0]        == (1/12)*L_5*(m.asin(L_1*m.sin(m.theta)/m.hs)+m.acos((m.hs**2+L_3**2-L_4**2)/(2*L_3*m.hs)))) # position of polished rod, inches 
    [m.Equation(m.v[i+1].dt()== m.a**2 * (m.u[i+2] - 2.0*m.u[i+1] + m.u[i])/m.dx**2 - m.pi*m.a*m.nu/(2.0*m.TVD)*m.v[i+1] - (1-m.rho_w*m.gamma/m.rho_r)*m.g) for i in range(npx-2) ]# wave equation
    m.Equation(m.q_out       == m.A_t * m.u[-1].dt()*12/231/42 * (1+m.y)/2) # rate of fluid production, barrels/
    
    # Equations for calculating rod loading
    # Load at surface
    m.Equation(m.f[0] == m.E*m.Ac*1/2/m.dx *(-m.u[2] + 4*m.u[1] -3*m.u[0]))
    # Load at pump
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

#MHE###########################################################################
#Parameters (Holds Measured values from MHE)  
fm = mhe.Param(value = sim.f[0])                                          

#SetGlobalOptions(MHE)##########################################################
mhe.options.IMODE = 5     # 4 = Dynamic Simulation (Seqential)
mhe.options.NODES = 2    # 3 = 3 Nodes, 2 = No collocation nodes
mhe.options.SOLVER = 3    # 1 =APOPT, 3 = IPOPT
mhe.options.time_shift = npt-1 # time shift forward for multiple simulations
mhe.options.MAX_ITER = 700
mhe.Obj((mhe.f[0] - fm)**2)

#SetLocalOptions (MHE)###############################################################
##FV #Variable to estimate
mhe.h.FSTATUS = 0
mhe.h.STATUS = 1
mhe.h.DMAX = 0.5

#MV
mhe.SPM_in.FSTATUS = 1
mhe.SPM_in.STATUS = 0

#Solve#########################################################################
#%%
# Solve the simulation in a loop to simulate a longer horizon
loops = 180 # number of steps forward in time

res = {}
solve_stat = np.zeros(loops)
t_cycle = 0
#PID Options ######################################################
e = np.zeros(loops)
ie = np.zeros(loops)
op = np.zeros(loops)

#SPM_output
op[0] = np.ones(1)*sim.SPM_in.value
ophi = 15
oplo = 5

pv = np.zeros(loops)
P = np.zeros(loops)
I = np.zeros(loops)
pv_ave = np.zeros(loops)

#Tuning
kc = 10
tauI = 5

###################################################################
#Initialize Storage Values
sim_ts = np.ones(npt)*sim.tsi.value # simulation time storage
sim_hstor = np.ones(npt)*sim.h.value # height of fluid in annulus storage
sim_q_ins= np.ones(npt)*sim.q_in.value # reservoir influx storage
sim_q_outs = np.ones(npt)*sim.q_out.value # production rate storage
sim_P_ress = np.ones(npt)*sim.P_res.value # reservoir pressure storage
sim_Vps    = np.ones(npt)*sim.Vp.value # cumulative volume produced storage
sim_NPVs   = np.ones(npt)*sim.NPV.value # NPV storage
sim_W_rods = np.ones(npt)*sim.W_rod.value # work of rod (work to lift fluid) storage
sim_SPMs   = np.ones(npt)*sim.SPM_in.value # Strokes per minute/ Torque storage Set Points
sim_SPMr   = np.ones(npt)*sim.SPM.value    #SPM storage
sim_thetas   = np.ones(npt)*sim.theta.value#Theta storage
sim_P_wfs  = np.ones(npt)*sim.P_wf.value # bottom hole pressure storage
sim_ys     = np.ones(npt)*sim.y.value # sign of du/dt storage
    

#MHE Storage
mpc_ts = np.empty(0) # simulation time storage
mhe_us = [np.array(mhe.u[i].value) for i in range(npx)] # u relative position storage
mhe_vs = [np.array(mhe.v[i].value) for i in range(npx)]
mhe_fs = [np.array(mhe.f[i].value) for i in range(npx)] # dynamic load storage 
mpc_hstor = np.empty(0)# height of fluid in annulus storage
mpc_q_ins= np.empty(0) # reservoir influx storage
mpc_q_outs = np.empty(0) # production rate storage
mpc_P_ress = np.empty(0) # reservoir pressure storage
mpc_Vps    = np.empty(0) # cumulative volume produced storage
mpc_NPVs   = np.empty(0) # NPV storage
mpc_W_rods =np.empty(0) # work of rod (work to lift fluid) storage
mpc_SPMs   = np.empty(0) # Strokes per minute/ Torque storage Set Points
mpc_SPMr   = np.empty(0)    #SPM storage
mpc_thetas   = np.empty(0)#Theta storage
mpc_P_wfs  = np.empty(0) # bottom hole pressure storage
mpc_ys     = np.empty(0) # sign of du/dt storage
mpc_Skins = np.empty(0) #Skin storage
mpc_porosity = np.empty(0)
mpc_A_d = np.empty(0)
mhe_Skins = np.empty(0) #Skin storage
mhe_hstor = np.empty(0)# height of fluid in annulus storage
mpc_Skins = np.empty(0)
mpc_Skinss = np.empty(0)
mpc_h = np.empty(0)
mpc_hss = np.empty(0)
###############################################################

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
        solve_stat[i] = t_cycle
    ##MHE##################################################################
    #Insert Measurements
    fm.value = sim.f[0].value
    
    #Insert move
    mhe.SPM_in.value = sim.SPM_in.value
    
    #Solve
    t_start = time.time()
    mhe.solve()

    #Pass values to MPC
    pv[i] = mhe.h.NEWVAL
    pv_ave[i] = np.average(pv[i-10:i])
    #Store new values for plotting
    mhe_hstor = np.append(mhe_hstor,mhe.h.value)
    
#    ##PID #####################################################################
#    if i<10:
    e[i] =  pv[i] - SP
#    if i>10:
#        e[i] = pv_ave[i] - SP
#    
#    if i >= 1:
    ie[i] = ie[i-1] + e[i]*tf
    
    P[i] = kc * e[i]
    I[i] = kc/tauI * ie[i]
     
    op[i] = op[0] + P[i] + I[i] 
    
    #Anti-reset Wind up
    if op[i] > 15:
        op[i] = ophi
        ie[i] = ie[i] - e[i]*tf
    if op[i] < 5:
        op[i] = oplo
        ie[i] = ie[i] - e[i]*tf
    
    #Pass output to simulation
    t_end   = time.time()
    t_cycle = t_end - t_start
    sim.SPM_in.value = np.ones(npt)*op[i]
    
#    #######################################################################
#    
    # Plotting
    plt.clf()
    ax=plt.subplot(311)
    ax.grid()
    plt.plot(sim_ts[0:i*npt],sim_SPMs[0:i*npt],'ro',label='SPM Set Point')
    plt.plot(sim_ts[0:i*npt],sim_SPMr[0:i*npt],'bo',label='SPM')
    plt.ylabel('Strokes per Minute')
    plt.legend(loc=2)
    ax=plt.subplot(312)
    ax.grid()
    plt.plot(sim_ts[0:i*npt],sim_hstor[0:i*npt],'k-',label= 'Height')
    plt.plot(sim_ts[0:i*npt], np.ones(i*npt)*SP, label = 'height SP')
    #plt.plot(ts[0:i*npt], mpc_hs[0:npt*i], label = 'mpc height')
    plt.ylabel('Annular Fluid Height (ft)')
    plt.legend(loc='best')
    ax = plt.subplot(313)
    ax.grid()
    plt.plot(sim_ts[0:i*npt], sim_q_outs[0:i*npt], label = 'q_out')
    plt.plot(sim_ts[0:i*npt], sim_q_ins[0:i*npt], label = 'q_in')
    plt.legend()
    plt.ylabel('Flow Rate, STB/s')
    plt.xlabel('Time (sec)')
    plt.draw()
    plt.pause(0.02)

#%%

res['solve_stat'] = solve_stat
res['ts'  ]     = sim_ts
res['us'  ]     = sim_us
res['vs'    ]   = sim_vs
res['fs']       = sim_fs
res['hstor'  ]  = sim_hstor
res[ 'q_ins' ]  = sim_q_ins
res[ 'q_outs' ] = sim_q_outs
res['P_ress'  ] = sim_P_ress
res['Vps'    ]  = sim_Vps
res[ 'NPVs'  ]  = sim_NPVs
res['W_rods' ]  = sim_W_rods   
res['SPMs']     = sim_SPMs
res['SPMr']     = sim_SPMr
res['thetas']   = sim_thetas
res['P_wfs']    = sim_P_wfs
res[ 'ys']      = sim_ys    
res['h_SP']     = np.ones(loops*npt)*SP

np.save('PI_MHE_Control_Results_Aggresive_20npt_10npt.npy', res) 
#%%
# Load dictionary of results
res = np.load('PI_MHE_Control_Results_Aggresive_20npt_10npt.npy').item()

#%% Plotting from dictionary


plt.figure()
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
ax=plt.subplot(311)
ax.grid()
plt.plot(res['ts'], res['SPMs'], 'r--', label=r'$T_{net}$ (ft-lb)')#, s = 4, c='b' ) # 'ro', for latex
plt.plot(res['ts'], res['SPMr'],'b-',  label=r'Actual')#, s = 4, c = 'r') #'bo',
plt.ylabel('SPM', fontsize = 12)
plt.legend(loc= 1,fontsize = 12)
plt.xlim(0,180)
ax=plt.subplot(312)
ax.grid()
plt.plot(res['ts'], res['hstor'],'k-',label= 'Actual')
plt.plot(res['ts'], np.ones(np.size(res['ts']))*(sim.TVD.value*3/4 -3), label = 'SP') # fix
#plt.plot(ts[0:i*npt], mpc_hs[0:npt*i], label = 'mpc height')
plt.ylabel('Fluid Level (ft)', fontsize = 12)
plt.legend(loc=1,fontsize = 12)
plt.xlim(0,180)
ax = plt.subplot(313)
ax.grid()
plt.plot(res['ts'], res['q_outs'], label = r'$q_{out}$')
plt.plot(res['ts'], res['q_ins'], label = r'$q_{in}$')
plt.legend(loc = 1,fontsize = 12)
plt.ylabel('Flow (STB/s)', fontsize = 12)
plt.xlabel('Time (seconds)', fontsize = 12)
plt.xlim(0,180)
plt.draw()
   
#plt.legend(fontsize = 12)
#plt.ylabel('Fluid Level (ft)', fontsize = 12)
#plt.xlabel('Time (seconds)', fontsize = 12)

plt.tight_layout()
plt.savefig('PI_MHE_Control_K_10_Tau_5.eps', transparent = True, dpi = 1200)
plt.show()
#%%
# timing figure
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.figure()
plt.plot(np.linspace(1,180,179), res['solve_stat'][1:],'r-',label='Solve Time')
plt.plot(np.linspace(1,loops,loops-1), np.ones(loops-1)*np.average(res['solve_stat'][1:]),'b--', label = 'Average')
plt.plot(np.linspace(1,loops,loops-1), np.ones(loops-1),'k:', label = 'Real Time')
plt.xlabel('Control Cycle', fontsize = 12)
plt.ylabel('Computation Time (seconds)', fontsize = 12)
plt.legend(fontsize = 12)
plt.ylim(0,4)
plt.xlim(0,180)
plt.savefig('PI_MHE_Simulation_Timing.eps', dpi = 1200, transparent = True)
plt.show()















#%%
## Figure 1 ## Height in the Annulus and incoming and outgoing flow rate####
plt.figure(1, figsize = (6,4.5))
plt.subplot(211)
plt.plot(sim_ts, sim_hstor, 'r--', label = 'height in annulus')
plt.plot(sim_ts, np.ones(len(sim_ts))*SP, 'b--', label = 'height Set Point')
plt.ylabel('height, ft')
plt.legend()

plt.subplot(212)
plt.plot(sim_ts, sim_q_ins, 'b--', label = r'$q_{in}$')
plt.plot(sim_ts, sim_q_outs, 'g--', label = r'$q_{out}$')
plt.ylabel('Flow Rate, STB/s')
plt.xlabel('time, sec')
plt.legend()

plt.show()

##Figure 2: Reservoir Pressure Decline and Cumulative Volume Produced####
plt.figure(2, figsize = (6,4.5))
plt.subplot(211)
plt.plot(sim_ts, sim_P_ress, 'k--', label = 'Reservoir Pressure')
#plt.plot(m.time, P_wf.value, 'r--', label = r'$P_{wf}$')
plt.ylabel('Pressure, psi')
plt.legend()

plt.subplot(212)
plt.plot(sim_ts, sim_Vps, '--', label = 'Cumulative Volume Produced')
plt.ylabel('Volume, STB')
plt.xlabel('time, sec')
plt.legend()
plt.tight_layout()

plt.show()

##Figure 3: NPV #########################################
plt.figure(3, figsize = (6,4.5))
plt.plot(sim_ts, sim_NPVs/(1e6), 'g:', label = 'NPV')
plt.xlabel('time, sec')
plt.ylabel('NPV, $ Millions')
plt.legend()

plt.show()

########################################################

#Figure 4# Work of Motor And 
plt.figure(4, figsize = (6,4.5))
plt.subplot(311)
plt.plot(sim_ts,sim_W_rods, 'b-', label = 'Work Supplied by Motor' )
plt.ylabel('KiloWatts, KW')

plt.subplot(312)
plt.plot(sim_ts, sim_SPMs, 'r-', label = 'Input' )
plt.ylabel('SPM')

plt.subplot(313)
plt.plot(sim_ts, sim_P_wfs, 'r--', label = r'$P_{wf}$')
plt.ylabel('FBHP, psi')
plt.xlabel('time, sec')
#plt.tight_layout()

plt.show()

##Figure 5 -Doublet Test
plt.figure(5, figsize = (6,4.5))
plt.subplot(211)
plt.plot(sim_ts, sim_hstor, 'r--', label = 'height in annulus')
plt.ylabel('height, ft')
plt.legend()

plt.subplot(212)
plt.plot(sim_ts, sim_SPMs, 'b--', label = r'SPM')
plt.ylabel('strokes/min')
plt.xlabel('time, sec')
plt.legend()

plt.show()

