% 2016 ESE650 Final Option 1 - Multiple Plane Detection
% data load example script
clear all;
close all;

D = imread('./terrain/depth/05.png'); % uint16 (mm)
rgb = imread('./terrain/rgb/05.png'); % uint8

% Remove out of range depths 
D(D>2800) = 0;
D(D<0) = 0;

% Median Filter to smooth noisy depths
% D = medfilt2(D,[3 3]);
% D = medfilt2(D,[4 4]);
% D = medfilt2(D,[5 5]);
D = medfilt2(D,[15 15]);

% Camera calibration
% intrinsic matrix K time 2D image coordinates to get world coordinates in 3D
% X = f(x/z), Y = f(y/z) where f is focal length
% provided by Bhoram Lee
[pcx, pcy, pcz, r, g ,b] = depthToCloud(D, rgb);

pts = [pcx pcy pcz];


%%

tic;

kd_tree = KDTreeSearcher(pts,'distance','euclidean');
k = 64;
n = knnsearch(kd_tree,pts,'k',(k+1));

n = n(:,2:end);
p = repmat(pts(:,1:3),k,1) - pts(n(:),1:3);
p = reshape(p, size(pts,1),k,3);

% Covariance
P = zeros(size(pts,1),6);
P(:,1) = sum(p(:,:,1).*p(:,:,1),2);
P(:,2) = sum(p(:,:,1).*p(:,:,2),2);
P(:,3) = sum(p(:,:,1).*p(:,:,3),2);
P(:,4) = sum(p(:,:,2).*p(:,:,2),2);
P(:,5) = sum(p(:,:,2).*p(:,:,3),2);
P(:,6) = sum(p(:,:,3).*p(:,:,3),2);
P = P ./ k;

normals = zeros(size(pts));

for i = 1:size(pts,1)
    A = [P(i,1) P(i,2) P(i,3); P(i,2) P(i,4) P(i,5); P(i,3) P(i,5) P(i,6)];
    
    [U,D,~] = svd(A);
    p = U(:,3)';
    normals(i,:) = p(1:3);
end

pts1 = pts - repmat([0,0,10],size(pts,1),1);
if 1
    [~,idx] = max(abs(normals),[],2);
    idx = (1:size(normals,1))' + (idx-1)*size(normals,1);
    dir = normals(idx).*pts1(idx) > 0;
else
    dir = sum(normals.*pts1,2) > 0;
end

normals(dir,:) = -normals(dir,:);
toc;

%%

tic;
%display_normals( rgb, pts, normals', r, g, b, 1.4 );
display_normals( rgb, pts, normals', r, g, b, 1.0);
toc

