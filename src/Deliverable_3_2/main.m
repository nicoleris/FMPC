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
x = [0; 0; 0; 0];
y = [0; 0; 0; 0];
z = [0; 0];
yaw = [0; 0];

x_ref = 1;
y_ref = 1;
z_ref = 1;
yaw_ref = pi/4;

% Get control inputs with
ux = mpc_x.get_u(x, x_ref);
uy = mpc_y.get_u(y, y_ref);
uyaw = mpc_yaw.get_u(yaw, yaw_ref);
uz = mpc_z.get_u(z, z_ref);

Tf = 40;
sim_x = simulation(Tf, Ts, sys_x, mpc_x, x, x_ref);
sim_y = simulation(Tf, Ts, sys_y, mpc_y, y, y_ref);
sim_z = simulation(Tf, Ts, sys_z, mpc_z, z, z_ref);
sim_yaw = simulation(Tf, Ts, sys_yaw, mpc_yaw, yaw, yaw_ref);

%%
sim = quad.sim(mpc_x, mpc_y, mpc_z, mpc_yaw);
quad.plot(sim);