function [sample_bvecs, sample_idx] = bd_subsample(bvecs,n)

%
% function [sample_bvecs, sample_idx] = bd_subsample(bvecs,n)
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

elec_points = bd_read_epoints(n);

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
sample_idx = zeros(1, n);

% We will knock these out as we go along:
potential_idx = [1:length(bvecs)];

% Loop over the rotated eq points: 
for i=1:n
    this = new_pts(:,i); 
    % For each one of the points, calculate the distance from all the bvecs: 
    delta = zeros(1,size(bvecs,2));  
    for j=1:length(bvecs)
        delta(j) = vector_angle(this, bvecs(:, j));
    end
    
    this_idx = find(delta==min(delta));
    
    % find and choose the bvec corresponding to the smallest distance from
    % where you ended up and add to the return value 
    sample_bvecs(:,i) = bvecs(:, this_idx);
    sample_idx(i) = potential_idx(this_idx); 

    % Knock out that index, so that it doesn't recur: 
    potential_idx = cat(2, potential_idx(1:this_idx-1),...
                           potential_idx(this_idx+1:end));

    bvecs = cat(2, bvecs(:, 1:this_idx-1), bvecs(:, this_idx+1:end));

end

