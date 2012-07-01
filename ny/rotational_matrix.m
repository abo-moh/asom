function [ R ] = rotational_matrix(theta)
% ROTATIONAL_MATRIX Generates matrix to rotate points/vectors
    R = [cosd(theta), sind(theta);-sind(theta) cosd(theta)];
end