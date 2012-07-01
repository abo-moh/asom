function [ g_mag_o, g_dir_o ] = gradients_filter( g_mag, g_dir, varargin )
%GRADIENTS_FILTER Filters the gradients so that gradients with appropriate
% directions are present
% Set default values for optional arguments
angle_min = pi;
angle_max = 3/2*pi;

% Set values of given optional arguments
c_vargs = length(varargin);
for i = 1:2:c_vargs
    switch varargin{i}
        case 'min' %Set minimum angle
            angle_min = varargin{i+1};
        case 'max' %Set maximum angle
            angle_max = varargin{i+1};
        otherwise
            %Invalid input!
    end
end
% Only gradients with angles between angle_min and angle_max are returned
% Create copies of gradient matrices
g_mag_o = g_mag;
g_dir_o = g_dir;

% Generate filter mask
g_mask = (g_dir >= angle_min) & (g_dir <= angle_max);

% Apply filter mask
g_mag_o(~g_mask) = false;
g_dir_o(~g_mask) = false;
end