function [ clusterCtrs, pts2Cluster, cluster2Pts ] = meanshift( data_pts, pts, sigma )

[~,numPts] = size(data_pts);

stopThresh = sigma*1e-1;
sigmaSq = sigma^2;
remaining_inds = 1:numPts;
remainingNumPts = numPts;

clusterCtrs = [];
visited = false(1,numPts);   

num_clusters = 0;
% alpha = 0.000001;
% alpha = 5e-4;
% alpha = 5e-5;
% alpha = 8e-4;
% alpha = 4e-6;
alpha = 1e-5;
%alpha = 6e-6;
beta = 5;

lastNum = 0;

while remainingNumPts
    
    lastNum = remainingNumPts;
    rand_ind = randperm(remainingNumPts,1);
    start_ind = remaining_inds(rand_ind);
    eucl_mean = pts(start_ind,1:3)';
    clust_mean = data_pts(:,start_ind);
    clust_members = [];  
    this_clust_votes = zeros(1,numPts,'uint16');
    
    fprintf('remaining points: %i\n', remainingNumPts);
        
    iter = 0;
    while iter < 300
        
        iter = iter+1;
        
        meanEuclPts = repmat(eucl_mean', numPts,1);
        meanNormPts = repmat(clust_mean,1,numPts);
        
        d = sum((meanEuclPts' - pts(:,1:3)').^2);
        %angle = 5*real(acos(dot(meanNormPts,data_pts)));
        angle = beta*real(acos(dot(meanNormPts,data_pts)));
        
        combinedDist = angle + alpha*d;
        inds = find(combinedDist < sigmaSq);
        
        if isempty(inds) 
            inds = start_ind;
            iter = 300;
            fprintf('empty neighbors\n');
        end
        
        % innovation
        % rotating points to xy-plane
        
        tangent0 = cross(clust_mean, [1;0;0]);
        if dot(tangent0,tangent0) < 0.001
            tangent0 = cross(clust_mean, [0;1;0]);
        end
        tangent0 = tangent0 / norm(tangent0);
        tangent1 = cross(clust_mean,tangent0);
        tangent1 = tangent1 / norm(tangent1);
        R = [tangent0 tangent1 clust_mean];
        
        X0 = (R'*(pts(inds,1:3)' - meanEuclPts(inds,:)'));
        ind3 = X0(3,:)>.5 | X0(3,:)<-.5;
        X0(3,ind3) = 0;
        
        %figure(1); clf; plot3(X0(1,:),X0(2,:),X0(3,:), 'bo');
        %drawnow;
        
        X0 = R*X0 + meanEuclPts(inds,:)';
        
%         figure(1); clf; plot3(X0(1,:),X0(2,:),X0(3,:), 'bo');
%         hold on;
%         uvw = eucl_mean+clust_mean;
%         quiver3(eucl_mean(1),eucl_mean(2),eucl_mean(3), uvw(1),uvw(2),uvw(3),1);
%         hold off;
%         drawnow;
%         axis equal;
%         pause;
        
        this_clust_votes(inds) = this_clust_votes(inds)+1;
        
        prior_mean = [clust_mean; alpha*eucl_mean];
        t = mean(atan2(data_pts(2,inds), data_pts(1,inds)));
        p = mean(acos(data_pts(3,inds)));
        x = sin(p) * cos(t);
        y = sin(p) * sin(t);
        z = cos(p);
        clust_mean = [x;y;z];
        clust_mean = clust_mean ./ norm(clust_mean);
        %eucl_mean = mean(pts(inds,1:3),1)';
        eucl_mean = mean(X0,2);                 % innovation
        
        clust_members = [clust_members inds];
        visited(clust_members) = 1;
        
        dst = norm([clust_mean; alpha*eucl_mean] - prior_mean);
        
        if dst < stopThresh    
            
            closestCluster = 0;
            for i = 1:num_clusters
                dist = beta*real(acos(dot(clust_mean, clusterCtrs(1:3,i)))) + alpha*sum((eucl_mean - clusterCtrs(4:6,i)).^2);
                %adist = acosd(dot(clust_mean, clusterCtrs(1:3,i)));
                %edist = alpha*sum((eucl_mean - clusterCtrs(4:6,i)).^2)
                if dist < sigma/2
                    closestCluster = i;
                    break;
                end
            end
            
            if closestCluster > 0
                t = atan2(clust_mean(2), clust_mean(1));
                p = acos(clust_mean(3));
                cl_clust = clusterCtrs(:,closestCluster);
                t2 = atan2(cl_clust(2), cl_clust(1));
                p2 = acos(cl_clust(3));
                t = (t+t2)/2;
                p = (p+p2)/2;
                x = sin(p) * cos(t);
                y = sin(p) * sin(t);
                z = cos(p);
                avg_mean = 0.5*(eucl_mean+clusterCtrs(4:6,closestCluster));
                clusterCtrs(:,closestCluster) = [x;y;z; avg_mean];
                clust_votes(closestCluster,:) = clust_votes(closestCluster,:) + this_clust_votes;
            else
                num_clusters = num_clusters+1;
                clusterCtrs(:,num_clusters) = [clust_mean; eucl_mean];
                clust_votes(num_clusters,:) = this_clust_votes;
            end
            break;
        end
    end
    
    remaining_inds = find(visited == 0);
    remainingNumPts = length(remaining_inds);
end

[~,pts2Cluster] = max(clust_votes,[],1);

cluster2Pts = cell(num_clusters,1);
for i = 1:num_clusters
    clust_members = find(pts2Cluster == i);
    cluster2Pts{i} = clust_members;
end


end

