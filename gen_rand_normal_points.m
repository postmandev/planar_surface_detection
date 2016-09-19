
for i = 1:10
    % Create a normal vector
    if i==1
        disp('Testing yz noisy planar points')
        test_normal = [1 0 0]';
    elseif i==2
        disp('Testing xz noisy planar points')
        test_normal = [0 1 0]';
    elseif i==3
        disp('Testing xy noisy planar points')
        test_normal = [0 0 1]';
    else
        disp('Testing random skew noisy planar points')
        test_normal = rand(3,1);
    end

    test_normal = test_normal / norm(test_normal);
    [x,y] = gen_unit_vectors(test_normal);

    % Generate random points in the plane
    pts = [];
    for j = 1:360
        pt = x * (rand() * 2 - 1);
        pt = pt + y * (rand() * 2 - 1);
        pt = pt + test_normal * (rand() * 0.2 - 0.1);
        pts(j,:) = pt';
    end
    test_result = compute_best_plane(pts, true);
    
    pause;
end

