function [J] = pendulum_Jacobian(c,t,u,v)

% Function to calculate Jacobian of generalized forces for ALE-ANCF pendulum model
%
% INPUTS:
% c - structure containing model parameters:
%       n - number of elements
%       axial_stiffness - axial cable stiffness EA (N)
%       transverse_stiffness - transverse/bending stiffness EI (Nm^2)
% t - time value (s)
% u - vector of generalized coordinates
% v - vector of generalized velocities
%
% OUTPUTS:
% J - Jacobian matrix
%

J = [];

for i=1:c.n
   qe = [u(7*(i-1)+1:7*(i-1)+6),u(7*(i-1)+8:7*(i-1)+13)];
   p1 = u(7*(i-1)+7);
   p2 = u(7*(i-1)+14);
   
   Ji = sparse(-ALEANCF_forceJacobian_axial(c.axial_stiffness,p1,p2,qe) - ...
        ALEANCF_forceJacobian_transverse(c.transverse_stiffness,p1,p2,qe));
   J = blkdiag(J,Ji); 
   
end

B = ALEANCF_embed(c.n); %Embedding matrix
J = B'*J*B;

end