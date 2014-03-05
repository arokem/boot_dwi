function elec_points = bd_read_epoints(n)

bd_dir = fileparts(which(mfilename));

pts_dir = fullfile(bd_dir,'camino_pts'); 

elec_points = dlmread(fullfile(pts_dir, sprintf('Elec%03d.txt', n))); 

% The first line should be equal to n: 
assert(elec_points(1)==n, 'There is something wrong with the camino points file'); 

% The format is: n,x1,y1,z1,x2,y2,z2 ... xn,yn,zn 
elec_points = [elec_points(2:3:end),... 
               elec_points(3:3:end),...    
               elec_points(4:3:end)];