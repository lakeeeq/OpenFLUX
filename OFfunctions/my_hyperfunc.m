function y = my_hyperfunc(x,par)

% if x < 0.99999999
%     y = x * par / (1 - abs(x));
% else
%     y = 99999999;
% end

% if x < 0.99999999
    y = x * par ./ (1 - abs(x));
% else
%     y = ones(size(x))*99999999;
% end