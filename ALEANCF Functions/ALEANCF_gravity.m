function Q = ALEANCF_gravityHHT(mass_per_unit_length,p1,p2,qe)

% Function to calculate generalized gravity force for an ALE-ANCF element
%
% INPUTS:
% mass_per_unit_length - cable linear density (kg/m)
% p1 - material coordinate of first element node
% p2 - material coordinate of second element node
% qe - vector of reduced generalized coordinates
%
% OUTPUTS:
% Q - generalized force vector
%

Q = mass_per_unit_length*...
    [0;
     0;
     (981*p1)/200 - (981*p2)/200;
     0;
     0;
     -(327*(p1 - p2)^2)/400;
     (981*qe(9))/200 - (981*qe(3))/200 - (327*p1*qe(6))/200 + (327*p2*qe(6))/200 + (327*p1*qe(12))/200 - (327*p2*qe(12))/200;
     0;
     0;
     (981*p1)/200 - (981*p2)/200;
     0;
     0;
     (327*(p1 - p2)^2)/400;
     (981*qe(9))/200 - (981*qe(3))/200 + (327*p1*qe(6))/200 - (327*p2*qe(6))/200 - (327*p1*qe(12))/200 + (327*p2*qe(12))/200];

end