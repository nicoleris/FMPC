clc
clear all
close all

x = [0; 0; 0; 2];
y = [0; 0; 0; 0];
z = [0; 0];
yaw = [0; 0];

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

% Get control inputs with
ux = mpc_x.get_u(x);
uy = mpc_y.get_u(y);
uz = mpc_z.get_u(z);
uyaw = mpc_yaw.get_u(yaw);
