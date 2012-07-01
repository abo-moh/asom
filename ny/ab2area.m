function [ i_bw ] = ab2area( width, height, slope, offset )
%AB2POINTS Summary of this function goes here
%   Detailed explanation goes here
    i_bw = poly2mask([0 0 (-1 * (offset+1) ./ slope)], [0 (offset+1) 0], width, height);
end