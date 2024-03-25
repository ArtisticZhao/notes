function []=plot_H(h_grid, postion, titlename)
% show_h 显示信道估计的结果
% h_grid
if nargin < 2
postion = 1:size(h_grid, 2);
end
if nargin < 3
    titlename= 'channel estimate';
else
    titlename = sprintf("channel estimate: %s", titlename);
end
nsym = length(postion);
t=tiledlayout (2, nsym);
layout = 1:2*nsym;
layout = reshape(layout, [], 2).';
title(t, titlename)
ax_abs=[];
ax_ang=[];
for i=1:nsym
    % 绘制幅度
    axx = nexttile(layout (1, i));
    plot (abs (h_grid (:,postion(i))))
    title(sprintf("mag of sym %d", postion(i)));
    ax_abs = [ax_abs axx];
    % 绘制角度
    axx = nexttile(layout(2, i));
    plot (rad2deg(angle(h_grid(:,postion(i)))), 'o')
    title(sprintf("ang of sym %d", postion(i)));
    ax_ang = [ax_ang axx];
end
linkaxes (ax_abs, 'xy');
linkaxes (ax_ang, 'xy');
end