function norm_x = unit_vector(x)

% function unit_vector(x)
%
% Normalized the input vector x, such that it lies on the unit sphere

norm_x = x/l2_norm(x);
