function [mx, my, c] = fn_calc_plane_through_point_cloud(x, y, z, force_plane_through_origin)
%calculates best fit plane through point cloud given by (x,y,z) coordinates
%and returns answer in terms of coefficients of equation z = mx * x + my *
%y + c. If force_plane_through_origin = 1 then c is forced to be zero.
if force_plane_through_origin
    M = [sum(x .^ 2), sum(x .* y); ...
        sum(x .* y), sum(y .^ 2)];
    V = [sum(z .* x); sum(z .* y)];
    tmp = M \ V;
    mx = tmp(1);
    my = tmp(2);
    c = 0;
else
    M = [sum(x .^ 2), sum(x .* y), sum(x); ...
        sum(x .* y), sum(y .^ 2), sum(y); ...
        sum(x), sum(y), length(z)];
    V = [sum(z .* x); sum(z .* y); sum(z)];
    tmp = M \ V;
    mx = tmp(1);
    my = tmp(2);
    c = tmp(3);
end
end