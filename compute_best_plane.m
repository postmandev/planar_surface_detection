function [ n,mu ] = compute_best_plane( pts, show_plot )
%COMPUTE_BEST_PLANE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    show_plot = false;
end

x = pts(:,1);
y = pts(:,2);
z = pts(:,3);
N = length(pts);

mu = [mean(x), mean(y), mean(z)];
offset = [x-mu(1), y-mu(2), z-mu(3)];
A = (1/N)*(offset'*offset);

[U,D,~] = svd(A);
s = D(3,3);
p = U(:,3)';
p(4) = -p*mu';

n = p(1:3);

if show_plot
    range = max(max(pts) - min(pts)) / 2;
    mid_pt = (max(pts) - min(pts)) / 2 + min(pts);
    xlim = [-1 1]*range + mid_pt(1);
    ylim = [-1 1]*range + mid_pt(2);
    zlim = [-1 1]*range + mid_pt(3);

    h = plot3(pts(:,1),pts(:,2),pts(:,3), 'bo');
    hold on;
    set(get(h,'Parent'), 'DataAspectRatio', [1 1 1], 'XLim',xlim,'YLim',ylim,'ZLim',zlim);

    norm_pts = [mean(pts); mean(pts) + (n*range)];
    h = plot3(norm_pts(:,1),norm_pts(:,2),norm_pts(:,3), 'r-', 'LineWidth',3);
    title(sprintf('Normal vector: <%0.3f, %0.3f, %0.3f>', n), 'FontWeight', 'bold', 'FontSize',14);
    grid on;
    axis equal;
    hold off;
end


end

