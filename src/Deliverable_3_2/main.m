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

%%
x = [0; 0; 0; 0];
x_ref = [0; 0; 0; 2];

% Get control inputs with
ux = mpc_x.get_u(x, x_ref);


% sim = quad.sim(mpc_x, mpc_y, mpc_z, mpc_yaw);
% quad.plot(sim);