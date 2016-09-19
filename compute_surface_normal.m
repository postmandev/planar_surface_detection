function [ n ] = compute_surface_normal( pts )
%COMPUTE_SURFACE_NORMAL Summary of this function goes here
%   Detailed explanation goes here

solutions = zeros(3,1);

for i = 1:3
    
    P = pts;
    P(:,i) = 1;
    P_dot = P' * P;
    if det(P_dot) == 0
        solutions(i) = 0;
        continue;
    end
    solutions(i) = 1;
    coeff = P_dot \ (P' * pts(:,i));
    coeff_neg = -coeff;
    coeff_neg(i) = 1;
    coeff(i) = 1;
    n(:,i) = coeff_neg / norm(coeff);
end

if sum(solutions) == 0
    return;
end

mu = mean(pts);
offset = [pts(:,1)-mu(1) pts(:,2)-mu(2) pts(:,3)-mu(3)];
for i = 1:3
    if solutions(i) == 0
        res_sum(i) = NaN;
        continue;
    end
    residuals = offset * n(:,i);
    res_sum(i) = sum(residuals .* residuals);
end

best_fit = find(res_sum == min(res_sum));

n = n(:, best_fit(1));

range = max(max(pts) - min(pts)) / 2;
mid_pt = (max(pts) - min(pts)) / 2 + min(pts);
xlim = [-1 1]*range + mid_pt(1);
ylim = [-1 1]*range + mid_pt(2);
zlim = [-1 1]*range + mid_pt(3);

h = plot3(pts(:,1),pts(:,2),pts(:,3), 'bo');
hold on;
set(get(h,'Parent'), 'DataAspectRatio', [1 1 1], 'XLim',xlim,'YLim',ylim,'ZLim',zlim);

norm_pts = [mean(pts); mean(pts) + (n'*range)];
h = plot3(norm_pts(:,1),norm_pts(:,2),norm_pts(:,3), 'r-', 'LineWidth',3);
title(sprintf('Normal vector: <%0.3f, %0.3f, %0.3f>', n), 'FontWeight', 'bold', 'FontSize',14);
grid on;
axis equal;
hold off;


end

