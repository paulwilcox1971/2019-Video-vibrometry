function [dx, dy] = fn_calc_image_shift_vs_time(mov, options)
default_options.debug_axes = [];
default_options.db_thresh = -40;
default_options.force_to_dominant_image_feature_direction = 0;
default_options.filter = 'none';
default_options.return_velocity = 0;
options = fn_set_default_fields(options, default_options);

n = size(mov);

%vectors and matrices for image-space and k-space
if n(2) > 1
    x = linspace(-1, 1, n(2));
    kx = [0:n(2) - 1] / n(2);
    kx_shift = [-floor(n(2) / 2):ceil(n(2) / 2) - 1] / n(2);
else
    x = 0;
    kx = 0;
    kx_shift = 0;
end
if n(1) > 1
    y = linspace(-1, 1, n(1));
    ky = [0:n(1) - 1] / n(1);
    ky_shift = [-floor(n(1) / 2):ceil(n(1) / 2)-1] / n(1);
else
    y = 0;
    ky = 0;
    ky_shift = 0;
end
[X, Y] = meshgrid(x, y);
[KX, KY] = meshgrid(kx, ky);
[KX_SHIFT, KY_SHIFT] = meshgrid(kx_shift, ky_shift);

R = sqrt(X .^ 2 + Y .^ 2);
KR = sqrt(KX .^ 2 + KY .^ 2);
KR_SHIFT = sqrt(KX_SHIFT .^ 2 + KY_SHIFT .^ 2);

%set up image filter
switch options.filter
    case 'none'
        filt = ones(n(1), n(2));
    case 'gaussian'
        filt = exp(-(R * sqrt(-log(10 ^ (options.db_thresh / 20)))).^ 2);
    case 'disc'
        if n(1) > 1 & n(2) > 1
            filt = double(R <= 1);
        else
            filt = ones(n(1), n(2));
        end
    case 'taper'
        if n(1) > 1 & n(2) > 1
            tmp = R;
        elseif n(1) > 1
            tmp = y';
        elseif n(2) > 1
            tmp = x;
        end
        taper_size = 0.5;
        filt = (1+cos((tmp-taper_size)/(1-taper_size)*pi)) / 2;
        filt(tmp <= taper_size) = 1;
        filt(tmp >= 1) = 0;
end

%2D FFTs of each image frame
kmov = zeros(size(mov));
for t = 1:n(3)
    mov(:,:,t) = mov(:,:,t) .* filt;
    kmov(:,:,t) = fftshift(fft2(mov(:,:,t)));
end

if options.force_to_dominant_image_feature_direction & n(1) > 1 & n(2) > 1
    tmp = mean(abs(kmov),3) .* exp(-(R * sqrt(-log(10 ^ (options.db_thresh / 20))).^ 2));
    cov_matrix = [tmp(:)' * KX_SHIFT(:) .^ 2, tmp(:)' * (KX_SHIFT(:) .* KY_SHIFT(:)); tmp(:)' * (KX_SHIFT(:) .* KY_SHIFT(:)), tmp(:)' * KY_SHIFT(:) .^ 2];
    [a,b] = eig(cov_matrix);
    theta = atan2(a(2,1), a(1,1)) + pi/2;
end

vx = zeros(1,n(3));
vy = zeros(1,n(3));
for t = 2:n(3)
    %calculate phase shift between consecutive k-images
    dphi = angle(kmov(:,:,t) ./ kmov(:,:,t-1));
    %check where both k-space images are >  threshold
    k_image_over_thresh = (abs(kmov(:,:,t)) / max(max(abs(kmov(:,:,t)))) > 10 ^ (options.db_thresh / 20)) .* ...
        (abs(kmov(:,:,t-1)) / max(max(abs(kmov(:,:,t-1)))) > 10 ^ (options.db_thresh / 20)) .* ...
        (KR_SHIFT < 0.25);
    %find max KR for which all abs(dphi(k_image_over_thresh)) < pi / 2;
    k_image_dphi_over_thresh = abs(dphi) > pi / 4;%should be an option
    max_KR = min(KR_SHIFT(k_image_over_thresh & k_image_dphi_over_thresh));
    if isempty(max_KR)
        max_KR = 0.25;
    end
    valid_pts = (KR_SHIFT <= max_KR) & k_image_over_thresh;
    valid_pt_indices = find(valid_pts);
    if n(1) > 1 & n(2) > 1
        [mx, my, c] = fn_calc_plane_through_point_cloud(KX_SHIFT(valid_pt_indices), KY_SHIFT(valid_pt_indices), dphi(valid_pt_indices), 1);
        if options.force_to_dominant_image_feature_direction
            tmp = [1, tan(theta); 1, -1/tan(theta)] \ [mx + my * tan(theta); 0];%would be better without tans!!
            mxx = tmp(1);
            myy = tmp(2);
            tmp = -tmp / (2*pi);
        else
            tmp = -[mx; my] / (2 * pi);
        end
        vx(t) = tmp(1);
        vy(t) = tmp(2);
        if ~isempty(options.debug_axes) && t == floor(n(3)/2) %second bit just to only show this at one time step
            tmp = squeeze(abs(mov(:,:,t)));
            tmp = 20*log10(tmp / max(tmp(:)));
            cla(options.debug_axes(1));
            imagesc(x, y, tmp, 'Parent', options.debug_axes(1));caxis(options.debug_axes(1), [-40, 0]);
            axis(options.debug_axes(1), 'xy');
            axis(options.debug_axes(1), 'equal');
            
            tmp = squeeze(abs(kmov(:,:,t)));
            tmp = 20*log10(tmp / max(tmp(:)));
            cla(options.debug_axes(2));
            imagesc(kx_shift, ky_shift, tmp, 'Parent', options.debug_axes(2));caxis(options.debug_axes(2), [-40, 0]);
            hold(options.debug_axes(2), 'on');
            viscircles(options.debug_axes(2), [0,0], max_KR);
            if options.force_to_dominant_image_feature_direction
                plot(options.debug_axes(2), [-cos(theta), cos(theta)], [-sin(theta), sin(theta)], 'r')
            end
            axis(options.debug_axes(2), 'xy');
            axis(options.debug_axes(2), 'equal');
            
            cla(options.debug_axes(3));
            plot3(options.debug_axes(3), KX_SHIFT(valid_pt_indices), KY_SHIFT(valid_pt_indices), dphi(valid_pt_indices), 'r.');
            hold(options.debug_axes(3), 'on');
            a = axis(options.debug_axes(3));
            xx = [a(1:2); a(1:2)]; yy = [a(3:4); a(3:4)]'; zz = mx*xx+my*yy; surf(options.debug_axes(3), xx,yy,zz);
            if options.force_to_dominant_image_feature_direction
                xx = [a(1:2); a(1:2)]; yy = [a(3:4); a(3:4)]'; zz = mxx*xx+myy*yy; surf(options.debug_axes(3), xx,yy,zz);
            end
            axis(options.debug_axes(3), 'equal');
        end
    elseif n(2) > 1
        [m, c] = fn_calc_line_through_point_cloud(KX_SHIFT(valid_pt_indices), dphi(valid_pt_indices), 1);
        vx(t) = -m / (2 * pi);
        vy(t) = 0;
    elseif n(1) > 1
        [m, c] = fn_calc_line_through_point_cloud(KY_SHIFT(valid_pt_indices), dphi(valid_pt_indices), 1);
        vx(t) = 0;
        vy(t) = -m / (2 * pi);
%         if t == floor(n(3)/2)
%             keyboard
%         end
    end
end
if options.return_velocity
    dx = vx;
    dy = vy;
else
    dx = cumsum(vx);
    dy = cumsum(vy);
end
end