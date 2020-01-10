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
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% YOUR CODE HERE %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%epsilon_speed = opti.variable(1,N+1);
f = @(x,u) quad.f(x, u);
f_discrete = @(x,u) RK4(x,u,h,f);

thrust_parameter=0.6;

%---- objective ---------
opti.minimize(...
  1*(X(10,:) - REF(1))*(X(10,:) - REF(1))' + ... % Go to ref X
  1*(X(11,:) - REF(2))*(X(11,:) - REF(2))' + ... % Go to ref Y
  20*(X(12,:) - REF(3))*(X(12,:) - REF(3))' + ... % Go to ref Z
  0.5*(X(6,:) - REF(4))*(X(6,:) - REF(4))'  + ... % Go to ref Yaw
  thrust_parameter*U(1,:)*U(1,:)' + ... % Minimize cmd
  thrust_parameter*U(2,:)*U(2,:)' + ... % Minimize cmd
  thrust_parameter*U(3,:)*U(3,:)' + ... % Minimize cmd
  thrust_parameter*U(4,:)*U(4,:)'); % Minimize cmd



for k=1:N % loop over control intervals
  opti.subject_to(X(:,k+1) == f_discrete(X(:,k), U(:,k)));
end

opti.subject_to(0 <= U <= 1.5);
opti.subject_to(X(:,1)==X0);


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