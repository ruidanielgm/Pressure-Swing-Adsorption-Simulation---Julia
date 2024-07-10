using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets;

#System Configurations:
yinfeed = [0.25, 0.75];
y_start = [1, 0];
t_max = 12000;

#System Parameters -> Some of these will just be used later when I add the energy balance
R=8.3144; # (dm3 . kPa) / (K . mol)
nocomp=2; # CH4, CO2
Mwi=[1.60e-2,4.40e-2,4.0e-3]; # Kg/mol
Phigh=1e5; #Pa
#Plow=1e5;
Tinlet=40+273.15; #K
Tinf=40+273.15; #K
    # Bed properties
L=0.236; #m
Dw=0.021; #m
ebed=0.37; 
row=8238; #Kg/m3
Cpw=500;
alfaw=430.66;
alfawl=451.95; 
rob=565.2; #Kg/m3
bed_area=3.1415*(Dw/2)^2; # PI.r^2 #m2
    #Feed conditions
q_feed_slpm= 0.2; #SLPM
q_feed = 0.2E-3; #m3/min Stand.
# Adsorbent properties 
Rp  =1e-3; #m
apM =3/Rp; #1/m
rop =rob/(1-ebed);   #Kg/m3
ep  = 0.5; 
dp  =2*Rp; #m
Cps =900; 
    # Adsorption equilibrium
w=[1.0,1.0]; #Equilibrium Tuning Parameters [-] #novo
    #ISO QUAD
qsat  =  [    7.15,   6.86]; #mol/kg
S0    =  [ 7.35e-4,2.08e-4] /1e5; #bar-1/1e5 = Pa-1
DH    =  [   12590,  18842]; #J/mol
    # Transport parameters
kldf     =[1e-0, 1e-0];    #0.0005*90, 0.00001*900
Dax       =4.42E-4; #m2/s
lambda    =2.10; # heat axial dispersion
kf        =2.33e-2; #m/s
hf        =364; 
U         =49; 
hw        =3240; 
visc      =1.56e-5;
    #Properties
a0 = [    4.568,    3.259, 2.500];
a1 = [ -8.975E-3, 1.36E-3,     0];
a2 = [  3.631E-5,  1.50E-5,     0];
a3 = [ -3.407E-8, -2.37E-8,     0];
a4 = [ 1.091E-11, 1.06E-11,     0];

T = Tinlet; #later T will vary when I add the energy balance. For now, let's consider its constant

@parameters z t

@variables y1(..) y2(..) Cg1(..) Cg2(..) Cgt(..) Cs1(..) Cs2(..) q1(..) q2(..) q1star(..) q2star(..) Pmp1(..) Pmp2(..) S1(..) S2(..) u0(..) P(..) Mw(..)

Dt = Differential(t)
Dz = Differential(z)
Dzz = Differential(z)^2

#System Boundaries
z_min = 0.0
t_min = 0.0
z_max = L
t_max = t_max;

#Domains
domains = [z ∈ Interval(z_min, z_max), t ∈ Interval(t_min, t_max)]


#Equations

eq =[
#Mass balance
    #Equation 1 - #Mass balance to comp i - PDE
    Dt(Cg1(z,t)) ~ (ebed*Dax*(Dzz(Cg1(z,t))) 
        - (Dz(u0(z,t))*Cg1(z,t) + u0(z,t)*Dz(Cg1(z,t)))
        - (1-ebed)*apM*kf*(Cg1(z,t)-Cs1(z,t)))  / ebed, 
    Dt(Cg2(z,t)) ~ (ebed*Dax*(Dzz(Cg2(z,t))) 
        - (Dz(u0(z,t))*Cg2(z,t) + u0(z,t)*Dz(Cg2(z,t))) 
        - (1-ebed)*apM*kf*(Cg2(z,t)-Cs2(z,t)))  / ebed,   
    #Equation 2 
    Cg1(z,t) ~ y1(z, t) * Cgt(z,t),
    Cg2(z,t) ~ y2(z, t) * Cgt(z,t),
    #Equation 3 
    y1(z,t) ~ 1 - y2(z,t),
    #Equation 4 
    Dt(q1(z,t)) ~ (apM*kf*(Cg1(z,t)-Cs1(z,t))) / rop,
    Dt(q2(z,t)) ~ (apM*kf*(Cg2(z,t)-Cs2(z,t))) / rop,
    #Equation 5 
    Dt(q1(z,t)) ~ kldf[1] * (q1star(z,t) - q1(z,t)),
    Dt(q2(z,t)) ~ kldf[2] * (q2star(z,t) - q2(z,t)),
    #Equation 6 
    Pmp1(z,t) ~ Cs1(z,t) * (R*T), #later, T will be T(z,t)
    Pmp2(z,t) ~ Cs2(z,t) * (R*T), #later, T will be T(z,t)
    #Equation 7 
    S1(z,t) ~ S0[1] * exp(DH[1]/(R*T)), #later, T will be T(z,t)  
    S2(z,t) ~ S0[2] * exp(DH[2]/(R*T)), #later, T will be T(z,t)
    #Equation 8 - Isotherm
    q1star(z,t) ~ w[1] * qsat[1] * S1(z,t) * Pmp1(z,t) / (1 + S1(z,t)*Pmp1(z,t) + S2(z,t)*Pmp2(z,t)),
    q2star(z,t) ~ w[2] * qsat[2] * S2(z,t) * Pmp2(z,t) / (1 + S1(z,t)*Pmp1(z,t) + S2(z,t)*Pmp2(z,t)), 
#Momentum Balance
    Mw(z,t) ~ Mwi[1] * y1(z,t) + Mwi[2] * y2(z,t), #Just Auxiliary 
    Dz(P(z,t)) ~ (150*visc*(1-ebed)^2/(ebed^3*dp^2)*u0(z,t) + 1.75*(1-ebed)/(ebed^3*dp)*Cgt(z,t)*Mw(z,t)*u0(z,t)*abs(u0(z,t)))   /(-1),
    Cgt(z,t) ~ P(z,t) / (R * T)
]

#BCs and initial values
T_start =  40 + 273.15;
Cg1_start = y_start[1]*Phigh/(R*T_start);
Cg2_start = y_start[2]*Phigh/(R*T_start);
Cs1_start = Cg1_start;
Cs2_start = Cg2_start;
S1_start = S0[1] * exp(DH[1]/(R*T_start));
S2_start = S0[2] * exp(DH[2]/(R*T_start));
Pmp1_start = Cs1_start * (R*T_start);
Pmp2_start = Cs2_start * (R*T_start);
q1_start = w[1] * qsat[1] * S1_start * Pmp1_start / (1 + S1_start*Pmp1_start + S2_start*Pmp2_start);
q2_start = w[2] * qsat[2] * S2_start * Pmp2_start / (1 + S1_start*Pmp1_start + S2_start*Pmp2_start);
Cgt_start =  Phigh / (R*T_start);
u0inlet=(q_feed/bed_area)*(Tinlet/273.15) * (1/60) * (1e5/Phigh); #m/s
Cinlet = Phigh /(R*Tinlet); # inlet pressure/R/Tin  #mol/m3
bcs = [
    Dz(Cg1(0,t)) ~ (yinfeed[1]*u0inlet*Cinlet - Cg1(0,t)*u0(0,t)) / (-ebed*Dax), #BC
    Dz(Cg2(0,t)) ~ (yinfeed[2]*u0inlet*Cinlet - Cg2(0,t)*u0(0,t)) / (-ebed*Dax), #BC
    Dz(Cg1(z_max,t)) ~ 0,     #BC
    Dz(Cg2(z_max,t)) ~ 0,     #BC
    y1(z,0) ~ y_start[1], #IC
    y2(z,0) ~ y_start[2], #IC
    Pmp1(z,0) ~ Pmp1_start, #IC
    Pmp2(z,0) ~ Pmp2_start, #IC
    Cg1(z,0) ~ Cg1_start, #IC
    Cg2(z,0) ~ Cg2_start, #IC
    S1(z,0) ~ S1_start, #IC
    S2(z,0) ~ S2_start, #IC
    q1(z,0) ~ q1_start, #IC
    q2(z,0) ~ q2_start, #IC
    Dt(q1(z,0)) ~ 0,
    Dt(q2(z,0)) ~ 0,
    Cgt(z,0) ~  Cgt_start, #IC
    P(z_max,t) ~ Phigh, #BC
    u0(0,t) ~ u0inlet * Cinlet / Cgt(0,t), #BC
    u0(z,0) ~ u0inlet #<-------------- Initial Guess
] 

@named pdesys = PDESystem(eq,bcs,domains,[z,t],[y1(z,t),y2(z,t),Cg1(z,t),Cg2(z,t),Cgt(z,t),Cs1(z,t),Cs2(z,t),
        q1(z,t),q2(z,t),q1star(z,t),q2star(z,t), Pmp1(z,t), Pmp2(z,t), S1(z,t), S2(z,t), P(z,t), u0(z,t), Mw(z,t)])

#Discretization
N = 10
order = 2

discretization = MOLFiniteDifference([z=>N], 
                                      t,
                                      #advection_scheme = WENOScheme(), 
                                      advection_scheme = UpwindScheme(),
                                      approx_order = order, 
                                      #grid_align = <your grid type choice>,
                                      #should_transform = <Whether to automatically transform the PDESystem (see below)>
                                      #use_ODAE = <Whether to use ODAEProblem>
                                      )
prob = discretize(pdesys, discretization)

# Convert the PDE problem into an ODE problem
println("Discretization:")
@time prob = discretize(pdesys,discretization);

#solver
save_time = 10
println("Solver:")
@time sol = solve(prob, FBDF(autodiff = false), saveat = save_time, maxiters = 1e5, abstol = 1e-10, reltol = 1e-8)

#Solution and Plots
println("Results:")
discrete_z = sol[z]
discrete_t = sol[t]
y1_solution = sol[y1(z, t)]
y2_solution = sol[y2(z, t)]
P_solution = sol[P(z,t)]
Cg1_solution = sol[Cg1(z, t)]
Cg2_solution = sol[Cg2(z, t)]
Cgt_solution = sol[Cgt(z, t)]
u0_solution = sol[u0(z, t)]
F1 = u0_solution .* Cg1_solution / ebed
F2 = u0_solution .* Cg2_solution / ebed
q1_solution = sol[q1(z,t)]
q2_solution = sol[q2(z,t)]

using Plots

#histories
h_y = plot(discrete_t, y1_solution[end, :], title="y History", xlabel="t", ylabel="y", label="y1(z=L,t)")
plot!(discrete_t, y2_solution[end, :], label="y2(z=L,t)")
display(h_y)
h_F = plot(discrete_t, F1[end, :], title="F History", xlabel="t", ylabel="y", label="F1(z=L,t)")
plot!(discrete_t, F2[end, :], label="F2(z=L,t)")
display(h_F)
h_C = plot(discrete_t, Cg1_solution[end, :], title="Cg History", xlabel="t", ylabel="Cg", label="Cg1(z=L,t)")
plot!(discrete_t, Cg2_solution[end, :], label="Cg2(z=L,t)")
plot!(discrete_t, Cgt_solution[end, :], label="Cgt(z=L,t)")
display(h_C)
h_u0 = plot(discrete_t, u0_solution[end, :], title="u0 history", xlabel="t", ylabel="u0", label="u0(z=L,t)")
display(h_u0)
h_P = plot(discrete_t, P_solution[end, :], title="P history", xlabel="t", ylabel="u0", label="P(z=L,t)")
display(h_P)

desired_time = 500
time_input = floor(Int, desired_time / save_time)

P_y = plot(discrete_z, y1_solution[:, time_input], title="Profile Inside the Column", xlabel="z", ylabel="y", label="y1(z,t=defined)")
plot!(discrete_z, y2_solution[:, time_input], label="y2(z,t=0)")
display(P_y)
P_q = plot(discrete_z, q1_solution[:, time_input], title="q", xlabel="z", ylabel="q", label="q1(z,t=defined)", ylims=(0,1.2))
plot!(discrete_z, q2_solution[:, time_input], label="q2(z,t=defined)")
display(P_q)
P_p = plot(discrete_z, P_solution[:, time_input], title="P profile", xlabel="z", ylabel="P", label="P(z,t=defined)")
display(P_p)
P_u0 = plot(discrete_z, u0_solution[:, time_input], title="u0 profile", xlabel="z", ylabel="P", label="u0(z=defined)")
display(P_u0)


