function B = ALEANCF_embed(n)

% Function to calculate embedding matrix for ALE-ANCF cable model
% Refer to C. Westin, 'Modelling and simulation of marine cables with 
%   dynamic winch and sheave contact', p. 61-63
%
% INPUTS:
% n - number of elements
%
% OUTPUTS:
% B - embedding matrix
%

%%Define Independent and Dependent Coordinates
independent = 1:7;
for i=1:n
	independent = [independent (8:14)+14*(i-1)];
end

dependent = [];
for i=1:(n-1)
	dependent = [dependent (1:7)+14*(i)];
end    
        
%%Embedding Matrices
B1 = eye(length(independent));
B=zeros(14*n,length(independent));
for i=1:length(independent)
	B(independent(i),:) = B1(i,:);
end
for i=1:length(dependent)
	B(dependent(i),:) = B(dependent(i)-7,:);
end  

end