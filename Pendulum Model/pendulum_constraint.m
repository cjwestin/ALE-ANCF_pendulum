function [C_q] = pendulum_constraint(c,t,u)

% Function to calculate constraint matrix for ALE-ANCF pendulum model
%
% INPUTS:
% c - structure containing model parameters:
%       n - number of elements
% t - time value (s)
% u - vector of generalized coordinates
%
% OUTPUTS:
% C_q - constraint Jacobian matrix
%

C_q = zeros(c.n+4,7*c.n+7);

%Pin constraint
C_q(1:3,:) = [eye(3),zeros(3,7*c.n+4)];

%Constraints on material flow (p2)
for i = 1:c.n+1
	C_q(i+3,7*i) = 1;
end

end