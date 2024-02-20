function [u_out,v_out] = solver_SIHHT(t_span,u_0,v_0,force_fun,mass_fun,constraint_fun,Jacobian_fun,vJacobian_fun,options)

% Semi-implicit HHT ODE solver, integrates the ALE-ANCF model from
% t_span(1) to t_span(end)
%
% Refer to C. Westin and R.A. Irani, 'Efficient semi-implicit numerical 
% integration of ANCF and ALE-ANCF cable models with holonomic
% constraints', Computational Mechanics, 2023
%
% INPUTS:
% t_span - array of time values
% u_0 - vector of initial generalized coordinates
% v_0 - vector of initial generalized velocities
% force_fun - function handle, outputs column vector of generalized forces
% mass_fun - function handle, outputs mass matrix
% constraint_fun - function handle, outputs kinematic constraint matrix
% Jacobian_fun - function handle, outputs Jacobian matrix of the
%                generalized force vector Q with respect to the generalized
%                coordinates u
% vJacobian_fun - function handle, outputs Jacobian matrix of the
%                 generalized force vector Q with respect to the generalized
%                 velocities v
% options - structure containing solver options:
%            h - solver time-step (s)
%            alpha - numerical parameter for the HHT method, controls the
%                    amount of numerical damping
%            beta, gamma - numerical parameters for the HHT method
%
% OUTPUTS:
% u_out - array containing the generalized coordinates at each time-step
% v_out - array containing the generalized velocities at each time-step
%


%Initialize generalized coordinates and velocities
u = u_0;
v = v_0;

ndof = length(u_0);                                                         %Number of degrees of freedom

%Predict initial generalized accelerations
[C_q] = constraint_fun(0,u);                                                %Constraint Jacobian matrix
M = mass_fun(0,u,v);                                                        %Mass matrix
Q = force_fun(0,u,v);                                                       %Generalized forces

sol = [M,C_q';C_q,zeros(size(C_q,1))]\[Q;zeros(size(C_q,1),1)];             %Solve linear system of equations
a = sol(1:ndof);                                                            %Generalized accelerations
L = sol(ndof+1:end);                                                        %Lagrange multipliers

%Initialize output arrays
u_out = zeros(length(u_0),length(t_span));
v_out = zeros(length(v_0),length(t_span));

for i = 1:length(t_span)

    t = t_span(i);        
    u_out(:,i) = u;
    v_out(:,i) = v;

    %Calculate generalized forces
    Q = force_fun(t,u,v);

    %Calculate position Jacobian
    Ju = Jacobian_fun(t,u,v);

    %Calculate velocity Jacobian
    if ~isempty(vJacobian_fun)
        Jv = vJacobian_fun(sys,t,u,v);    
    else
        Jv = zeros(size(Ju));
    end
    
    %Calculate mass matrix
    M = mass_fun(t,u,v);

    %Calculate constraint Jacobian
    [C_q] = constraint_fun(t,u);

    %Form LHS matrix
    H = [1/(1+options.alpha)*M - options.gamma*options.h*Jv - ...
            options.beta*options.h^2*Ju,...
            C_q';C_q, zeros(size(C_q,1))];
    
    %Form RHS vector
    R = [options.alpha/(options.alpha+1)*(C_q'*L - Q) + Q + ...
        Ju*(options.h*v + 0.5*options.h^2*(1-2*options.beta)*a) + ...
        Jv*(1-options.gamma)*options.h*a;...
        zeros(size(C_q,1),1)];
    
    %Solve linear system of equations
    sol = H\R;
    
    a_new = sol(1:ndof);                                                    %Updated acceleration
    L_new = sol(ndof+1:end);                                                %Updated Lagrange multipliers

    %Update generalized coordinates and velocities
    u = u + options.h*v + 0.5*options.h^2*((1-2*options.beta)*a ...
                                                + 2*options.beta*a_new);
    v = v + (1-options.gamma)*options.h*a + options.gamma*options.h*a_new;
    
    %Values for next time-step
    a = a_new;
    L = L_new;

end

end