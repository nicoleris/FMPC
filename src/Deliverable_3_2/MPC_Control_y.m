classdef MPC_Control_y < MPC_Control
  
  methods
    % Design a YALMIP optimizer object that takes a steady-state state
    % and input (xs, us) and returns a control input
    function ctrl_opt = setup_controller(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   x(:,1) - initial state (estimate)
      %   xs, us - steady-state target
      % OUTPUTS
      %   u(:,1) - input to apply to the system
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      [n,m] = size(mpc.B);
      
      % Steady-state targets (Ignore this before Todo 3.2)
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % SET THE HORIZON HERE
      N = 300;
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 

      % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are 
      %       the DISCRETE-TIME MODEL of your system

      alpha_max = 0.035;
      Ma_max = 0.3;
      F = [1; -1]; f = [alpha_max; alpha_max];
      M = [1; -1]; m = [Ma_max; Ma_max];
      
      R = 1;
      Q = 1*eye(n);
      
      syst = LTISystem('A', mpc.A, 'B', mpc.B);
      syst.x.max = [inf; alpha_max; inf; inf];
      syst.x.min = [-inf; -alpha_max; -inf; -inf];
      syst.x.penalty = QuadFunction(Q);
      syst.u.penalty = QuadFunction(R);
      Qf = syst.LQRPenalty.weight;
      Yf = syst.LQRSet;
      Ff = Yf.A;
      ff = Yf.b;
      
      % WRITE THE CONSTRAINTS AND OBJECTIVE HERE
      con = [];
      obj = 0;
      
      
      con = con + (x(:,2) == mpc.A*(x(:,1) - xs) + mpc.B*(u(:,1) - us));
      con = con + (M*(u(1, 1) - us) <= m);
      obj = obj + (x(:, 1) - xs)'*Q*(x(:, 1) - xs) + (u(:,1) - us)'*R*(u(:,1) - us);
      
      for i = 2:N-1
          con = con + (x(:, i+1) == mpc.A*(x(:, i) - xs) + mpc.B*(u(:, i) - us));
          con = con + (F*(x(2, i) - xs(2, 1)) <= f) + (M*(u(1, i) - us) <= m);
          obj = obj + (x(:, i) - xs)'*Q*(x(:, i) - xs) + (u(:, i) - us)'*R*(u(:, i) - us);
      end
%       
%       con = con + (Ff*x(:,N) <= ff);
%       obj = obj + x(:,N)'*Qf*x(:,N);
%       
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      ctrl_opt = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
        {x(:,1), xs, us}, u(:,1));
    end
    
    
    % Design a YALMIP optimizer object that takes a position reference
    % and returns a feasible steady-state state and input (xs, us)
    function target_opt = setup_steady_state_target(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   ref    - reference to track
      % OUTPUTS
      %   xs, us - steady-state target
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Steady-state targets
      n = size(mpc.A,1);
      xs = sdpvar(n, 1);
      us = sdpvar;
      
      % Reference position (Ignore this before Todo 3.2)
      ref = sdpvar;            
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
      u_limit = 0.3;
      
      con = [-u_limit <= us <= u_limit, xs == mpc.A*xs + mpc.B*us, ref == mpc.C*xs];
      obj = us^2;
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
      
    end
  end
end
