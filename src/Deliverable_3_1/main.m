clc
clear all
close all

%% Setup
Ts = 1/5;
quad = Quad(Ts);
[xs, us] = quad.trim();
sys = quad.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_yaw] = quad.decompose(sys, xs, us);

% Design MPC controller
mpc_x = MPC_Control_x(sys_x, Ts);
mpc_y = MPC_Control_y(sys_y, Ts);
mpc_z = MPC_Control_z(sys_z, Ts);
mpc_yaw = MPC_Control_yaw(sys_yaw, Ts);

%% Simulate
x0 = [0; 0; 0; 2];
y0 = [0; 0; 0; 2];
z0= [0; 2];
yaw0 = [0; pi/4];

Tf = 8;
sim_x = simulation(Tf, Ts, sys_x, mpc_x, x0);
sim_y = simulation(Tf, Ts, sys_y, mpc_y, y0);
sim_z = simulation(Tf, Ts, sys_z, mpc_z, z0);
sim_yaw = simulation(Tf, Ts, sys_yaw, mpc_yaw, yaw0);
