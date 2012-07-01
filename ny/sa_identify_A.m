function [ A, ind ] = sa_identify_A( line_points )
%SA_IDENTIFY_A Identify A in vector of points

    % Identify potential points
    y_value = max(line_points(:,1));
    A_potentials = find(line_points(:,1) == y_value);
    
    % In case of more than one potential point, choose the middle point
    % (ceil in case of even number of potential points)
    ind = A_potentials(ceil(length(A_potentials)/2));
    A = line_points(ind,:);

end

