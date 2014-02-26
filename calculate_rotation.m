function rot = calculate_rotation(a,b)

% function rot = calculate_rotation(a,b)
%
% Calculate the rotation matrix from the unit vector a to 
% the unit vector b, up to sign inversions.  
%  
% Example: 
%
% >> a = [1 0 0]; 
% >> b = [0 1 0]; 
% >> rot = calculate_rotation(a,b); 
% >> rot * b'
% 
% ans =
% 
%     1.0000
%     0.0000
%          0
% 

% What is the angle of rotation between them? 
alpha = vector_angle(a,b); 

% For no rotation, return the identity:
if alpha==0
    rot = eye(3); 
    return 
end
% Find the orthonormal basis for the null-space  of these two 
% vectors. This is a unit vector orthogonal to both of them: 
u = null([a;b]); 

% Using the quaternion notation: 
% (See http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#From_the_rotations_to_the_quaternions)

w = cos(alpha/2); 
xyz = u * sin(alpha/2); 

% Now derive the rotation matrix: 
rot = quat2dcm([w xyz']); 


