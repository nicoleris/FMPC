function [ctrl, traj] = ctrl_NMPC(quad)
import casadi.*
opti = casadi.Opti(); % Optimization problem
N = 20; % MPC horizon [SET THIS VARIABLE]
h=0.2;
%decision variables
X = opti.variable(12,N+1); % state trajectory variables
U = opti.variable(4, N); % control trajectory (thrusts)
X0 = opti.parameter(12,1); % initial state
REF = opti.parameter(4,1); % reference position [x,y,z,yaw]

f = @(x,u) quad.f(x, u);
f_discrete = @(x,u) RK4(x,u,h,f);

param=[0.7 20 0.5 0.3];

%---- objective ---------
opti.minimize(...
  param(1)*(X(10,:)-REF(1))*(X(10,:)-REF(1))' + ... % tracking ref X
  param(1)*(X(11,:)-REF(2))*(X(11,:)-REF(2))' + ... % tracking ref Y
  param(2)*(X(12,:)-REF(3))*(X(12,:)-REF(3))' + ... % tracking ref Z
  param(3)*(X(6,:) - REF(4))*(X(6,:)-REF(4))' + ... % tracking ref Yaw
  param(4)*U(1,:)*U(1,:)' + ... % Minimize input
  param(4)*U(2,:)*U(2,:)' + ... % Minimize input
  param(4)*U(3,:)*U(3,:)' + ... % Minimize input
  param(4)*U(4,:)*U(4,:)');     % Minimize input

for k=1:N % loop over control intervals
  opti.subject_to(X(:,k+1) == f_discrete(X(:,k), U(:,k)));
end

opti.subject_to(0 <= U <= 1.5);
opti.subject_to(X(:,1)==X0);
opti.subject_to(X(10,:)<=3.5);
opti.subject_to(X(11,:)<=1);
opti.subject_to(X(12,:)<=REF(3));

ctrl = @(x,ref) eval_ctrl(x, ref, opti, X0, REF, X, U);
end

function u = eval_ctrl(x, ref, opti, X0, REF, X, U)
%Set the initial state and reference
opti.set_value(X0, x);
opti.set_value(REF, ref);
%Setup solver NLP
ops = struct('ipopt', struct('print_level',0, 'tol', 1e-3), 'print_time', false);
opti.solver('ipopt', ops);
%Solve the optimization problem
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');
u = opti.value(U(:,1));
% Use the current solution to speed up the next optimization
opti.set_initial(sol.value_variables());
opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end