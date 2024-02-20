function [M] = pendulum_mass(c,t,u,v)

% Function to calculate mass matrix for ALE-ANCF pendulum model
%
% INPUTS:
% c - structure containing model parameters:
%       n - number of elements
%       mass_per_unit_length - cable linear density (kg/m)
% t - time value (s)
% u - vector of generalized coordinates
% v - vector of generalized velocities
%
% OUTPUTS:
% M - ALE-ANCF mass matrix
%

M = [];

for i=1:c.n
   qe = [u(7*(i-1)+1:7*(i-1)+6),u(7*(i-1)+8:7*(i-1)+13)];
   p1 = u(7*(i-1)+7);
   p2 = u(7*(i-1)+14);

   M = blkdiag(M,ALEANCF_mass(c.mass_per_unit_length,p1,p2,qe)); 
end

%Embedding
B = ALEANCF_embed(c.n); %Embedding matrix
M = B'*M*B;

end