function sample_bvecs = bd_subsample(bvecs,n)

%
% function [x,y,z] = bd_sample(bvecs,n)
% 
% Generate a sub-sample of size n of directions from the provided x0,y0,z0. 
% 
% Parameters
% ----------
% bvecs: int array (3 by n), a set of cartesian coordinates for a set of
% bvecs 
% n: int, how many bvecs to sub-sample from this set. 
% 
% Returns 
% -------
% x,y,z: The coordinates of the sub-sample
% 
% Notes
% -----
% Directions are chosen from the camino-generated electro-static repulsion
% points in the directory camino_pts.
% 

bd_dir = fileparts(which(mfilename));

pts_dir = fullfile(bd_dir,'camino_pts'); 

elec_points = dlmread(fullfile(pts_dir, sprintf('Elec%03d.txt', n))); 

% The first line should be equal to n: 
assert(elec_points(1)==n, 'There is something wrong with the camino points file'); 

% The format is: n,x1,y1,z1,x2,y2,z2 ... xn,yn,zn 
elec_points = [elec_points(2:3:end),... 
               elec_points(3:3:end),...    
               elec_points(4:3:end)];

% Since the camino points cover only a hemi-sphere, a random half of the
% points need to be inverted by 180 degrees to account for potential
% gradient asymmetries:
idx = randperm(size(elec_points,1));
idx = idx(1:end/2); 
elec_points(idx,:) = elec_points(idx,:) * -1;  

% The seed bvec is the one relative to which the rest are chosen: 
seed = ceil(rand * length(bvecs));
seed_coords = bvecs(:,seed); 

pt1 = elec_points(1,:); 

% Find the rotation transformation from the first eq point to the seed: 
rot = calculate_rotation(seed_coords', pt1);

% The rotation is accurate up to a sign reversal in each coordinate, so we 
% check and correct:  
rot_back = rot * pt1';

% This has -1 for sign-change and 1 for no sign change: 
sign_changer = sign(-1 *((sign(rot_back) ~= sign(seed_coords)) - 0.5)); 

% Multiply each line separately: 
rot = [rot(1,:) * sign_changer(1);...
       rot(2,:) * sign_changer(2);...
       rot(3,:) * sign_changer(3)];
   
% Now apply the rotation to all the eq points: 
new_pts = rot * elec_points';

% Pre-allocate for the return values 
sample_bvecs = zeros(3,n); 

% Loop over the rotated eq points: 
for i=1:n
    this = new_pts(:,i); 
    % For each one of the points, calculate the distance from all the bvecs: 
    delta = zeros(1,size(bvecs,2));  
    for j=1:size(bvecs,2)
        delta(j) = vector_angle(this, bvecs(:,j));
    end
    % find and choose the bvec corresponding to the smallest distance from
    % where you ended up and add to the return value 
    sample_bvecs(:,i) = bvecs(:,delta==min(delta));
end





