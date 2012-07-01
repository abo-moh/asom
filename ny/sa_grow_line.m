function [ line_final, start_index ] = sa_grow_line( points, start )
%SA_GROW_LINE Filters away all points that aren't on the line defined by start index
%   The program starts at the "start" point, grows upwards until no more
%   points are added. Afterwards the program grows the line downwards from
%   the "start" point until no more points are added.

    line_points = [];
    threshold = 10;
    line_points{1} = start;
    line_points{2} = [];
    % Select points above or at same height as start point
    upper_points = points(...
                    points(:,2) <= start(2) & ...
                    (points(:,1) ~= start(1) | points(:,2) ~= start(2)),:);
    
    % Select points below start point
    lower_points = points(points(:,2) > start(2),:);
    
    grow_points = {upper_points, lower_points};
    for p = 1:length(grow_points)
        new = 1;
        current = start;
        while ~isempty(new)
            % Create mask to choose pixels within threshold perimeter
            mask = sqrt((grow_points{p}(:,1) - current(1)).^2 + ...
                        (grow_points{p}(:,2) - current(2)).^2) ...
                        <= threshold;
            new = find(mask);

            if ~isempty(new)
                % Find closest point to add to line if more than 1 possible
                % points
                if length(new) > 1
                    t_values = sqrt((grow_points{p}(new,1) - current(1)).^2 + ...
                                    (grow_points{p}(new,2) - current(2)).^2);
                    [~,min_ind] = min(t_values);
                    new = new(min_ind);
                end
                
                new_point = grow_points{p}(new,:);
                line_points{p}(end+1,:) = new_point;
                % Set new endpoint
                current = line_points{p}(end,:);
                % Set new point to 0 in order to
                grow_points{p}(new,1) = 0;
                grow_points{p}(new,2) = 0;
            end
        end
    end
    % Concatenate lines together into one line with upper startpoint
    % in line_final(1,:) and lower end point in line_final(end,:)
    line_final = [line_points{1}(end:-1:1,:);line_points{2}(1:end,:)];
    start_index = size(line_points{1},1);
    
    %% Debugging plots
    
    %{
    figure;% plot(upper_points(:,1), upper_points(:,2), 'Xr');
    hold on, axis equal, axis ij;
    plot(lower_points(:,1), lower_points(:,2), 'Xr');
    plot(line_final(:,1),line_final(:,2),'Xb');
    lp = line_final(end,:);
    up = line_final(1,:);
    cp = [lp; up];
    plot(start(1), start(2),'Db');
    plot(cp(1,1), cp(1,2), 'Pr');
    plot(cp(2,1), cp(2,2), 'Pg');
    %}
    
    
end

