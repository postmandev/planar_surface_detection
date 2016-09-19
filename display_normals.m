function display_normals( rgb, pts, n, r, g, b, sigma )

display1 = 0;

if display1
    figure
    pcshow(pts, [r' g' b']/255,'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down')
    title('Estimated normals of a point cloud')
    hold on

    x = pts(1:10:end, 1);
    y = pts(1:10:end, 2);
    z = pts(1:10:end, 3);
    u = n(1,1:10:end)';
    v = n(2,1:10:end)';
    w = n(3,1:10:end)';

    % Plot the normal vectors
    quiver3(x, y, z, u, v, w);
    hold off

    % Flip the normals to point towards the sensor location. This is only
    % necessary when the inward/outward direction of surface needs to be
    % determined.
    sensorCenter = [0, -0.3, 0.3]; % in X,Y,Z coordinates
    for k = 1 : numel(x)
       p1 = sensorCenter - [x(k), y(k), z(k)];
       p2 = [u(k), v(k), w(k)];
       % Flip the normal vector if it is not pointing towards the sensor
       angle = atan2(norm(cross(p1, p2)), p1*p2');
       if angle > pi/2 || angle < -pi/2
           u(k) = -u(k);
           v(k) = -v(k);
           w(k) = -w(k);
       end
    end
end


%tic
[clustCent,point2cluster,clustMembsCell] = meanshift(n, pts,sigma);
%toc

numClust = length(clustMembsCell);
clust_area = 400;

figure(10),clf,hold on
rand_colors = rand(numClust,3);
cnum = 0;
for k = 1:numClust
    myMembers = clustMembsCell{k};
    if numel(myMembers) < clust_area
        continue;
    end
    cnum = cnum+1;
    myClustCen = clustCent(:,k);
    plot3(n(1,myMembers),n(2,myMembers),n(3,myMembers), '.', 'Color', rand_colors(k,:))
    %plot3(myClustCen(1),myClustCen(2),myClustCen(3), 'o','MarkerEdgeColor','k','MarkerFaceColor',rand_colors(k,:), 'MarkerSize',10)
end
title(['Number of Clusters:' int2str(cnum)]); % int2str(numClust)
axis equal
drawnow;

figure

pcshow(pts, [r', g', b']/255,'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down')
%drawnow;
hold on;
for i = 1:numClust
    myMembers = clustMembsCell{i};
    if isempty(myMembers) || numel(myMembers) < clust_area
        continue;
    end
    pcshow(pts(myMembers,:), rand_colors(point2cluster(myMembers),:),'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down');
    drawnow;
end
%imshow(im);
hold off;
drawnow;

end

