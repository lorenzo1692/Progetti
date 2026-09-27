function fig = plot_conductor_zoom(out, t, fig_title)
%PLOT_CONDUCTOR_ZOOM Mesh and interfaces of one conductor of the FEM surrogate.
%
%   fig = PLOT_CONDUCTOR_ZOOM(out, t, fig_title) zooms the WP_MECH_SURROGATE
%   mesh on turn t: elements coloured by material (cable, jacket, turn
%   insulation, filler, inter-layer insulation), element edges including
%   the curved Q8/T6 midside nodes, the analytic cable and jacket outlines
%   of WP_TURN_GEOMETRY drawn over it, and the duplicated cable/jacket
%   contact nodes (red). For a RIS conductor the cable must appear as a
%   circle and the jacket thickness must grow from the side middles to
%   the corners.

if nargin < 2 || isempty(t), t = 1; end
if nargin < 3 || isempty(fig_title), fig_title = sprintf('Turn %d: mesh and interfaces', t); end
m = out.mesh; xy = m.xy; geo = out.geo;
cell = geo.cells(t,:);
pad = 0.15*(cell(2) - cell(1));
bb = [cell(1)-pad, cell(2)+pad, cell(3)-pad, cell(4)+pad];

fig = figure('Name', 'FEM conductor zoom', 'NumberTitle', 'off');
hold on
cols = [0.20 0.45 0.85; 0.55 0.55 0.55; 0.60 0.90 0.95; 0.85 0.85 0.85; 0.75 0.75 0.60; 0.85 0.35 0.15; 0.90 0.80 0.45];
qi = any(inbox(xy, m.q8, bb), 2);
ti = any(inbox(xy, m.t6, bb), 2);
for mm = 1:7
    e = find(qi & m.q8_mat(:) == mm);
    if ~isempty(e)
        patch('Faces', m.q8(e, [1 5 2 6 3 7 4 8]), 'Vertices', xy, 'FaceColor', cols(mm,:), ...
            'EdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.3);
    end
    e = find(ti & m.t6_mat(:) == mm);
    if ~isempty(e)
        patch('Faces', m.t6(e, [1 4 2 5 3 6]), 'Vertices', xy, 'FaceColor', cols(mm,:), ...
            'EdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.3);
    end
end
% analytic outlines
c = geo.cable_rr(t,:); j = geo.jacket_rr(t,:);
[x, y] = rr(c); plot(x, y, 'b-', 'LineWidth', 1.2);
[x, y] = rr(j); plot(x, y, 'k-', 'LineWidth', 1.2);
% contact pairs (jacket side node of each duplicated pair) of this turn
if isfield(m, 'kpair') && ~isempty(m.kpair)
    k = m.kpair(:,1);
    in = xy(k,1) > cell(1) & xy(k,1) < cell(2) & xy(k,2) > cell(3) & xy(k,2) < cell(4);
    plot(xy(k(in),1), xy(k(in),2), 'r.', 'MarkerSize', 10);
end
axis equal; axis(bb); box on
xlabel('Toroidal x [m]'); ylabel('Radial y [m]');
title({fig_title, sprintf('%s conductor: blue = analytic cable, black = jacket outer, red = cable/jacket contact nodes', ...
    geo.shape_name)}, 'Interpreter', 'none');
end

function tf = inbox(xy, E, b)
X = reshape(xy(E,1), size(E)); Y = reshape(xy(E,2), size(E));
tf = X > b(1) & X < b(2) & Y > b(3) & Y < b(4);
end

function [x, y] = rr(c)
% outline of [xa xb ya yb r]
xc = (c(1)+c(2))/2; yc = (c(3)+c(4))/2; w = c(2)-c(1); h = c(4)-c(3); r = c(5);
cx = [w/2-r, -w/2+r, -w/2+r, w/2-r]; cy = [h/2-r, h/2-r, -h/2+r, -h/2+r];
x = []; y = [];
for q = 1:4
    th = (q-1)*pi/2 + linspace(0, pi/2, 24);
    x = [x, xc + cx(q) + r*cos(th)]; %#ok<AGROW>
    y = [y, yc + cy(q) + r*sin(th)]; %#ok<AGROW>
end
x(end+1) = x(1); y(end+1) = y(1);
end
