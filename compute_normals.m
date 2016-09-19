function [normal,mu] = compute_normals(ptCloud,k)

X = ptCloud.Location;
[m,n] = size(X);
Dk = ones(n,1);
NN = ones(m,n);

normal = zeros(n,m);
mu = zeros(n,m);

for t = 1:1:m;
    x = X(t,:);
%     x = ones(m,1)*x;
%     dist = (X - x).^2;
%     dist = sqrt(sum(dist,2)');
%     dist = sort(dist);
%     Dk = dist(k+1);
%     %index = find(dist == Dk); 
%     NN = X(dist < Dk,:)
    [indices,dists] = findNearestNeighbors(ptCloud,x,k);
    [normal(:,t),mu(:,t)] = compute_best_plane(X(indices,:), false);
end

end
