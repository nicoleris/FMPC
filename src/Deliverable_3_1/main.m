clc
clear all
close all

%% Part 2
% quad = Quad();
% [xs,us] = quad.trim(); % Compute steady-state for which 0 = f(xs,us)
% sys = quad.linearize(xs, us); % Linearize the nonlinear model
% sys_transformed = sys * inv(quad.T); % New system is A * x + B * inv(T) * v
% [sys_x, sys_y, sys_z, sys_yaw] = quad.decompose(sys, xs, us);

%% Part 3
% discrete_system = c2d(sys, 1/5);
x = [0; 0; 0; 2];
y = [0; 0; 0; 2];
z = [0; -2];
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

sim = quad.sim(mpc_x, mpc_y, mpc_z, mpc_yaw);
quad.plot(sim);
