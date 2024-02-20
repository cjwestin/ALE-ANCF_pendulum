function [Q] = pendulum_forces(c,t,u,v)

% Function to calculate generalized forces for ALE-ANCF pendulum model
%
% INPUTS:
% c - structure containing model parameters:
%       n - number of elements
%       mass_per_unit_length - cable linear density (kg/m)
%       axial_stiffness - axial cable stiffness EA (N)
%       transverse_stiffness - transverse/bending stiffness EI (Nm^2)
% t - time value (s)
% u - vector of generalized coordinates
% v - vector of generalized velocities
%
% OUTPUTS:
% Q - column vector of generalized forces
%

Q = [];

for i=1:c.n
   qe = [u(7*(i-1)+1:7*(i-1)+6),u(7*(i-1)+8:7*(i-1)+13)];
   p1 = u(7*(i-1)+7);
   p2 = u(7*(i-1)+14);

   %qe_t = [v(7*(i-1)+1:7*(i-1)+6),v(7*(i-1)+8:7*(i-1)+13)];
   %p1_t = v(7*(i-1)+7);
   %p2_t = v(7*(i-1)+14);
   
   Q_gravity = ALEANCF_gravity(c.mass_per_unit_length,p1,p2,qe);
   Q_axial = ALEANCF_axialforce(c.axial_stiffness,p1,p2,qe);
   Q_transverse = ALEANCF_transverseforce(c.transverse_stiffness,p1,p2,qe);

   Q = [Q;Q_gravity - Q_axial - Q_transverse];
end

%Embedding
B = ALEANCF_embed(c.n); %Embedding matrix
Q = B'*(Q);


end