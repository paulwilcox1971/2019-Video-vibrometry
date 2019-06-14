function [m, c] = fn_calc_line_through_point_cloud(x, y, force_line_through_origin)
%calculates best fit line through point cloud given by (x,y) coordinates
%and returns answer in terms of coefficients of equation y = m * x + c.
%If force_line_through_origin = 1 then c is forced to be zero.
if force_line_through_origin
    m = sum(y .* x) / sum(x .^ 2);
    c = 0;
else
    M = [sum(x .^ 2), sum(x); ...
        sum(x), length(y)];
    V = [sum(y .* x); sum(y)];
    tmp = M \ V;
    m = tmp(1);
    c = tmp(2);
end
end