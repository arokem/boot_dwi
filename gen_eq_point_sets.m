% gen_eq_point_sets
% 
% Make and save eq point sets with http://eqsp.sourceforge.net/ for the 3d 
% sphere. This requires that the eqsp toolbox be available on the matlab
% path. 


for n=2:150 
    
    file_name = fullfile('eq_points', sprintf('eq_point_set%03d.txt', n)); 
    a = eq_point_set(2,n);   
    dlmwrite(file_name, a, 'newline','unix'); 

end
