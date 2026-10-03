function h = draw_turn_section(tg, k, x0, y0, col_tins, col_jacket, col_cable)
%DRAW_TURN_SECTION Draw one turn cell with its real conductor shape.
%
%   h = DRAW_TURN_SECTION(tg, k, x0, y0, col_tins, col_jacket, col_cable)
%   draws the turn of layer k whose cell has its lower-left corner at
%   (x0, y0), using the per-layer geometry tg from WP_TURN_GEOMETRY:
%   insulated cell, jacket (outer rounded rectangle, corner radius jk_r)
%   and cable (Rect: rounded rectangle; RIS: circle of diameter d). Returns
%   the three patch handles [insulation jacket cable].

w = tg.cell_w(k); hgt = tg.cell_h(k); ti = tg.tins;
xc = x0 + w/2; yc = y0 + hgt/2;
h(1) = patch([x0 x0+w x0+w x0], [y0 y0 y0+hgt y0+hgt], col_tins, 'EdgeColor', [0.3 0.3 0.3]);
[xj, yj] = rrect_outline(xc, yc, tg.jk_w(k), tg.jk_h(k), tg.jk_r(k));
h(2) = patch(xj, yj, col_jacket, 'EdgeColor', 'none');
[xs, ys] = rrect_outline(xc, yc, tg.cab_w(k), tg.cab_h(k), tg.cab_r(k));
h(3) = patch(xs, ys, col_cable, 'EdgeColor', 'none');
if w - 2*ti <= 0, set(h(2), 'Visible', 'off'); end
end

function [x, y] = rrect_outline(xc, yc, w, h, r)
% Closed outline of a w x h rectangle with corner radius r (r = w/2 = h/2
% gives a circle).
n = 16;
cx = [w/2-r, -w/2+r, -w/2+r, w/2-r]; cy = [h/2-r, h/2-r, -h/2+r, -h/2+r];
x = []; y = [];
for q = 1:4
    th = (q-1)*pi/2 + linspace(0, pi/2, n);
    x = [x, xc + cx(q) + r*cos(th)]; %#ok<AGROW>
    y = [y, yc + cy(q) + r*sin(th)]; %#ok<AGROW>
end
end
