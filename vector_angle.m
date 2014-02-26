function alpha = vector_angle(a,b)
% 
%function alpha  = vector_angle(a,b)
% 
% Find the angle between two vectors a and b


% Normalize them both to unit length: 
norm_a = unit_vector(a); 
norm_b = unit_vector(b); 

% If the vectors are the same, the angle is per-definition 0: 
if all(norm_a==norm_b)
    alpha = 0;  
else
    alpha = acos(dot(norm_a,norm_b)); 
end
