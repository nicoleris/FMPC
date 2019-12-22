clc
clear all
close all

Ts = 1/5;
quad = Quad(Ts);
[xs, us] = quad.trim();
sys = quad.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_yaw] = quad.decompose(sys, xs, us);

% Design MPC controller
mpc_x = MPC_Control_x(sys_x, Ts);
mpc_y = MPC_Control_y(sys_y,Ts);
mpc_yaw = MPC_Control_yaw(sys_yaw,Ts);
mpc_z = MPC_Control_z(sys_z,Ts);


%%  

x = [0; 0; 0; 2];
y = [0; 0; 0; 2];
z = [0; 2];
yaw = [0; pi/4];

Tf = 20;
sim_x = simulation(Tf, Ts, sys_x, mpc_x, x);
sim_y = simulation(Tf, Ts, sys_y, mpc_y, y);
sim_z = simulation(Tf, Ts, sys_z, mpc_z, z);
sim_yaw = simulation(Tf, Ts, sys_yaw, mpc_yaw, yaw);
