% Rescaled Bliss Independence
% Following method by Yeh's group
% Tekin et al. Ecology Letters. 2020
% Using a newly introduced framework to measure ecological stressor interactions

function [DA_xy, cat] = RBI(w_x, w_y, w_xy)

%w_x = 0.70; % effect of stressor 1
%w_y = 0.50; % effect of stressor 2
w_xy_pred = w_x * w_y; % expected effect

%w_xy = 0.6; % actual amount

% Case of synergy
if w_xy - w_xy_pred <= 0 % negative sign
    DA_xy = (w_xy - w_xy_pred) ./ abs(0 - w_xy_pred);

elseif w_xy - w_xy_pred > 0

    % Case of antagonistic buffering
    if w_xy < min(w_x, w_y)
        DA_xy = (w_xy - w_xy_pred) ./ abs(min(w_x, w_y) - w_xy_pred);
          
    % Case of antagonistic suppression
    elseif w_xy > min(w_x, w_y)
        DA_xy = 1 + (w_xy - min(w_x, w_y)) ./ abs(1 - min(w_x, w_y));
    end

end

% Print result
D_crit = 0.5; % threshold for criticality
if abs(DA_xy) < D_crit
    cat = 'additivity';
elseif DA_xy < 0
    cat = 'synergy';
elseif DA_xy > 0 & DA_xy <= 1
    cat = 'antagonistic buffering';
elseif DA_xy > 0 & DA_xy > 1
    cat = 'antagonistic suppression';
else
    cat = 'Unknown'
end

end % end of function