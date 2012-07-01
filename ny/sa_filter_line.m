function [ filtered_line_points ] = sa_filter_line( line_points, start )
%SA_FILTER_LINE Filter out points that don't belong to the line
%   Filters out linepoints added above the leftmost point on the linepart
%   over the start point. And filters out linepoints added below the leftmost
%   point on the linepart under the right-most point on the breast-boundary
%   below the startpoint

    %% Find upper bound x-value
    val_upper = min(line_points(1:start,1));
    % In case of more than one lowest x-values, choose the top one
    % (lowest y-value)
    potentials = find(line_points(:,1) == val_upper);
    [~, ind] = min(line_points(potentials,2));
    bound_upper = potentials(ind);
    
    %% Find lower bound x-value
    % Start by moving down the line until we've reached the max x-value
    % we do this by checking if any of the 10 neighbour pixels (5 on each
    % side of the point) have a higher x-value than the current one
    % x-value
    
    % Initialize variables
    delta_points = 10;
    current = start+ceil(delta_points/2);
    % Indexes of neighbour pixels
    neighbours = [(current - ceil(delta_points/2)),(current + ceil(delta_points/2))];
    
    % Repeat until max x-value of neighbour pixels is found
    while line_points(current,1) < max(line_points(neighbours(1):neighbours(2),1))
        % Find index of neighbour with higher x-value
        [~, next_max] = max(line_points(neighbours(1):neighbours(2),1));
        % Update to current to the neighbour just found and set new
        % neighbour indexes
        current = neighbours(1)-1+next_max;
        neighbours = [(current - ceil(delta_points/2)),(current + ceil(delta_points/2))];
    end
    
    val_lower = min(line_points(current:end,1));
    
    % In case of more than one lowest x-values, choose the bottom one
    % (highest y-value)
    potentials = find(line_points(:,1) == val_lower);
    [~, ind] = max(line_points(potentials,2));
    bound_lower = potentials(ind);
    
    %% Set filtered return linepoints
    filtered_line_points = line_points(bound_upper:bound_lower,:);
    
    %% Debugging plots
    
    %{
    figure, plot(line_points(:,1),line_points(:,2), 'Xr');
    hold on, axis equal, axis ij;
    plot(filtered_line_points(:,1), filtered_line_points(:,2),'Xb');
    plot(line_points(bound_upper,1),line_points(bound_upper,2),'Xg');
    plot(line_points(bound_lower,1),line_points(bound_lower,2),'Xg');
    %figure, plot(filtered_line_points(:,1),filtered_line_points(:,2),'Xb');
    %}
    
end

