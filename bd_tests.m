% Provides a test-suite for the functions used here:

% vector_angle (also tests other functions, such as unit_vector, l2_norm,
% etc).

% Orthogonal vectors:
a = [1 0 0];
b = [0 1 0];

assertEqual(vector_angle(a,b), pi/2);

% Equal vectors:
a = [randn, randn, randn];
b = a;
% Should be an angle of 0:
assertEqual(vector_angle(a,b), 0);

% Same direction:
a = [1 0 0];
b = [0.5 0 0];
assertEqual(vector_angle(a, b), 0);

% Opposite direction:

a = [1 0 0];
b = [-0.5 0 0];

assertEqual(vector_angle(a,b), pi);

% a = unit_vector([randn randn randn]);
% b = unit_vector([randn randn randn]);

a = [1 0 0];
b = [1 0 0];
rot = calculate_rotation(a,b);

assertEqual(a', rot *b')

a = [1 0 0];
b = [0 1 0];
rot = calculate_rotation(a,b);

assertAlmostEqual(a', rot *b')

a = [1 0 0];
b = [0 0 1];
rot = calculate_rotation(a,b);

assertAlmostEqual(a', rot *b')


% Test that the bvectors you get back are unique:
for i=1:10
    bv_orig = bd_read_epoints(100);
    bv_orig = bv_orig';
    [xyz, idx] = bd_subsample(bv_orig, 30);
    assertEqual(length(unique(idx)), length(idx));
end

assertEqual(bv_orig(:, idx), xyz);


