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
      N = 20;
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 

      % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are 
      %       the DISCRETE-TIME MODEL of your system

      %CONSTRAINTS PARAMETERS
      alpha_max = 0.035;
      Ma_max = 0.3;
      
      M = [1; -1]; m = [Ma_max; Ma_max];
      F = [0 1 0 0 ; 0 -1 0 0 ]; f = [alpha_max; alpha_max];
      
      R = 15;
      Q = diag([10,10,10,10]);      
      
      [K, Qf, ~] = dlqr(mpc.A, mpc.B, Q, R);
      K = -K;
      
      Yf = polytope([F; M*K], [f; m]);
      Acl = [mpc.A + mpc.B*K];
      while 1
        prevYf = Yf;
        [T,t] = double(Yf);
        preXf = polytope(T*Acl,t);
        Yf = intersect(Yf, preXf);
        if isequal(prevYf, Yf)
            break
        end
      end
      [Ff,ff] = double(Yf);
      
      
      %CONSTRAINTS AND OBJECTIVE
      con = [];
      obj = 0;

      con = con + (x(:,2) == mpc.A*x(:,1) + mpc.B*u(:,1));
      con = con + (M*u(:, 1) <= m);
      obj = obj + x(:, 1)'*Q*x(:, 1) + u(:,1)'*R*u(:,1);
      
      for i = 2:N-1
          con = con + (x(:, i+1) == mpc.A*x(:, i) + mpc.B*u(:, i));
          con = con + (F*x(:, i) <= f) + (M*u(:, i) <= m);
          obj = obj + x(:, i)'*Q*x(:, i) + u(:, i)'*R*u(:, i);
      end
      
      con = con + (Ff*(x(:,N) - xs) <= ff);
      obj = obj + (x(:,N) - xs)'*Qf*(x(:,N) - xs);
      
      
      %PLOT OF INVARIANT SET
      figure;
      
      subplot(3, 1, 1);
      Yf.projection(1:2).plot();
      xlabel('$\dot{\beta}$ [rad/s]','interpreter','latex')
      ylabel('$\beta$ [rad]','interpreter','latex');
      title('Dimensions 1 and 2');
      
      subplot(3, 1, 2);
      Yf.projection(2:3).plot();
      xlabel('$\beta$ [rad]','interpreter','latex')
      ylabel('$\dot{y}$ [m/s]','interpreter','latex');
      title('Dimensions 2 and 3');
      
      subplot(3, 1, 3);
      Yf.projection(3:4).plot();
      xlabel('$\dot{y}$ [m/s]','interpreter','latex')
      ylabel('$y$ [m]','interpreter','latex');
      title('Dimensions 3 and 4');
      
      sgtitle('Terminal invariant set for y');
      
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
      
      con = [];
      obj = 0;
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
      
    end
  end
end
