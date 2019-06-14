restoredefaultpath;
addpath(genpath('C:\Users\mepdw\Git\NDT-Library')); %change this to location of NDT library on your machine
clear;
close all;

%File details, image and time ranges to analyse

original_fname = 'tuning fork 3.avi'; 
%If rescan_video == 1, Matlab will open original_fname, crop, downsample and 
%save as Matlab file with same name but *.mat extension instead. 
%If rescan_video == 0, Matlab will open the *.mat file with same name and 
%use pre-cropped and downsampled data that has been saved.
rescan_video = 1; 

%Following only meaningful if rescan_video = 1
x_range = [380, 1280];y_range = [1, 300];
frame_range = [1, inf];
spatial_downsample = 1;

%Processing options
options.filter = 'gaussian';
options.force_to_dominant_image_feature_direction = 1;
options.return_velocity = 1;
sub_image_size_x = 41;
sub_image_size_y = 41;
xstep = max([floor(sub_image_size_x / 1),1]);
ystep = max([floor(sub_image_size_y / 1),1]);

show_debug_axes = 1;

%--------------------------------------------------------------------------
%load video file - either rescan raw file and save in Matlab form or use
%pre-saved Matlab one to save re-scanning.

[dummy, fname, dummy] = fileparts(original_fname);
mat_fname = [fname, '.mat'];
if rescan_video
    vidObj = VideoReader(original_fname);
    stored_frame_index = 0;
    read_frame_index = 0;
    while hasFrame(vidObj)
        tmp = sum(readFrame(vidObj), 3);
        if isempty(y_range)
            y_range = [1, size(tmp, 1)];
        end
        if isempty(x_range)
            x_range = [1, size(tmp, 2)];
        end
        read_frame_index = read_frame_index + 1;
        if read_frame_index >= frame_range(1) & read_frame_index <= frame_range(2)
            stored_frame_index = stored_frame_index + 1;
            m(stored_frame_index).data = tmp(y_range(1): spatial_downsample: y_range(2), x_range(1): spatial_downsample: x_range(2));
        end
    end
    %3D cube of downsampled data
    time_step = 1 / vidObj.FrameRate;
    mov = reshape([m(:).data], size(m(1).data, 1), size(m(1).data, 2), []);
    save(mat_fname, 'mov', 'time_step');
else
    load(mat_fname);
    xsize = size(mov,2);
    ysize = size(mov,1);
    tsize = size(mov,3);
end

%--------------------------------------------------------------------------
%show first frame of video so you can see what's what
figure;
imagesc(mov(:,:,1));colormap gray; axis equal;

figure;
if show_debug_axes
    options.debug_axes(1) = subplot(3,3,3);
    options.debug_axes(2) = subplot(3,3,6);
    options.debug_axes(3) = subplot(3,3,9);
end

axes1 = subplot(1,3,1);
axes2 = subplot(1,3,2);
imagesc(mean(mov, 3), 'Parent', axes1);
colormap(axes1, gray);
caxis(axes1, [0,1024]);
axis(axes1, 'equal');
axis(axes1, 'xy');
hold(axes1, 'on');
%--------------------------------------------------------------------------
h = plot(axes1, [0,0], [1,1], 'r.-');

xc = floor(sub_image_size_x / 2)+1:xstep:xsize-floor(sub_image_size_x / 2);
yc = floor(sub_image_size_y / 2)+1:ystep:ysize-floor(sub_image_size_y / 2);

dx_1d = zeros(length(xc), length(yc), tsize);
dy_1d = zeros(length(xc), length(yc), tsize);
dx_2d = zeros(length(xc), length(yc), tsize);
dy_2d = zeros(length(xc), length(yc), tsize);

%--------------------------------------------------------------------------
%the main loop - working through sub-images in video and extracting
%time-histories of each sub-image
for xi = 1:length(xc)
    xi_start = xc(xi) - floor(sub_image_size_x / 2);
    xi_end = xc(xi) + floor(sub_image_size_x / 2);
    for yi = 1:length(yc)
        yi_start = yc(yi) - floor(sub_image_size_y / 2);
        yi_end = yc(yi) + floor(sub_image_size_y / 2);
        delete(h);
        h = plot(axes1, [xi_start,xi_start,xi_end,xi_end,xi_start], [yi_start,yi_end,yi_end,yi_start,yi_start], 'r.-');
        
        %1D version
%         [dx_1d(xi, yi,:), dummy] = fn_calc_image_shift_vs_time(mov(yc(yi):yc(yi), xi_start:xi_end, :), options);
%         [dummy, dy_1d(xi, yi,:)] = fn_calc_image_shift_vs_time(mov(yi_start:yi_end, xc(xi):xc(xi), :), options);
        %2D version
        [dx_2d(xi, yi,:), dy_2d(xi, yi,:)] = fn_calc_image_shift_vs_time(mov(yi_start:yi_end, xi_start:xi_end, :), options);
        
        cla(axes2);
        plot(axes2, squeeze(dx_2d(xi, yi,:)), 'r.-');
        hold(axes2, 'on');
        plot(axes2, squeeze(dy_2d(xi, yi,:)), 'g.-');
        drawnow;
    end
end

%--------------------------------------------------------------------------
%various bits of post-processing to look at frequency content
dx_2d_spec = fft(dx_2d, 2048, 3) * 2 / size(dx_2d, 3);
dy_2d_spec = fft(dy_2d, 2048, 3) * 2 / size(dy_2d, 3);
mean_spec= squeeze(mean(mean(abs(dx_2d_spec),1),2))+squeeze(mean(mean(abs(dy_2d_spec),1),2));
[dummy,index_of_peak] = max(mean_spec(1:length(mean_spec)/2));
mean_phase = angle(mean(mean(dx_2d_spec(:,:,index_of_peak))) + mean(mean(dy_2d_spec(:,:,index_of_peak))));
x2d = real(dx_2d_spec(:,:,index_of_peak) * exp(-1i * mean_phase))';
y2d = real(dy_2d_spec(:,:,index_of_peak) * exp(-1i * mean_phase))';

figure;
imagesc(mov(:,:,1));
axis equal; 
axis tight; 
colormap gray;
hold on; 
[mxc, myc] = meshgrid(xc, yc); 
h = quiver(mxc, myc, x2d, y2d);
set(h, 'Linewidth', 2, 'Color', 'r');

figure;
subplot(2,1,1);
imagesc(mxc(:,1), myc(1,:), y2d);
colorbar;
title({sprintf('Subimage: %i x %i', sub_image_size_x, sub_image_size_y), ...
    sprintf('Steps: %i, %i', xstep, ystep)});
subplot(2,1,2);
plot(mxc(yi,:)', y2d');

