clear;clc;
addpath('ALEANCF Functions');
addpath('Pendulum Model');

%% Model Parameters
%Physical properties
c.diameter              = 0.01;                                             %(m)
c.area                  = pi/4*c.diameter^2;                                %(m^2)
c.moment_of_inertia     = 1e-8;                                             %(m^4)
c.elastic_modulus       = 2e7;                                              %(Pa)
c.density               = 5000;                                             %(kg/m^3)

c.mass_per_unit_length  = c.area*c.density;                                 %(kg/m)
c.axial_stiffness       = c.elastic_modulus*c.area;                         %(N)
c.transverse_stiffness  = c.moment_of_inertia*c.elastic_modulus;            %(Nm^2) 

%Mesh properties
c.n                     = 10;                                               %Number of elements
c.len_nom               = 0.1;                                              %Nominal element length, unstretched (m)

load('initial_coordinates.mat')                                             %Get initial generalized coordinate vector

%% Simulation options

%Numerical parameters for SI-HHT Solver
options.alpha = -0.05;                                                      %Controls amount of numerical damping      
options.beta = 0.25*(1-options.alpha)^2;
options.gamma = 0.5 - options.alpha;

options.h = 1e-3;                                                           %Time step (s)

t = 0:options.h:10;                                                         %Time vector (s)

%Function definitions
force_fun       = @(t,u,v)  pendulum_forces(c,t,u,v);
mass_fun        = @(t,u,v)  pendulum_mass(c,t,u,v);
constraint_fun  = @(t,u)    pendulum_constraint(c,t,u);
Jacobian_fun    = @(t,u,v)  pendulum_Jacobian(c,t,u,v);

%% Run simulation

tic;
[u,v] = solver_SIHHT(t,u_0,v_0,force_fun,mass_fun,constraint_fun,Jacobian_fun,[],options);
comptime = toc;

