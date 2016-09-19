function [ x,y ] = gen_unit_vectors( test_normal )

x = zeros(3,1);
while sum(x) == 0
    x = cross(test_normal,rand(3,1));
end
x = x / norm(x);
y = cross(test_normal,x);
y = y / norm(y);

if abs(dot(test_normal,x)) > 1e-8 || abs(dot(test_normal,y)) > 1e-8 || ...
        abs(dot(x,y)) > 1e-8
    error('There seems to be a problem ...');
end

end

