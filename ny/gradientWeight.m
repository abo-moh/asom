function [ weight ] = gradientWeight( gradients, a_hist, x, y )
%W Weight function used by Hough transformation
    % Returns weight calculated as
    % W(g(x,y)) = 1 / a_hist * integral(gradients_hist, 0,g(x,y))
    mult = 1/a_hist;
    integral = zeros(numel(x),1);
    for i = 1:numel(x)
        integral_t = histc(gradients(:),[0 gradients(y(i),x(i))]);
        integral(i) = integral_t(1);
    end
    weight = mult * integral;
end