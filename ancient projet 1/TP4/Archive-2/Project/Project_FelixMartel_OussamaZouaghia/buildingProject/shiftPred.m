%INPUT:
%  t - simulation time-step
%  N - Prediction Horizon of your MPC controller

%OUTPUTS:
% pred - prediction of the disturbance input, over the MPC prediction horizon
% cp - shifted price of electricity consumption over the MPC prediction horizon, at time step t
% sb - shifted comfort constraint off-sets over the MPC prediction horizon, at time step t


function [pred, cp, sb] = shiftPred(t, N, refCost, refSetbacks)

load building.mat;


%% Disturbance Prediction
pred = refDist(:,t:t+N-1);


%% Variable Price Prediction 
% *) Setting the costRef

cp = refCost(:,t:t+N-1);

%replace in section 3

%% Night-Setback Prediction

sb = refSetbacks(:,t:t+N-1);

%replace in section 4


end

