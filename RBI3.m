% DA = deviation from additivity
function [DA_R, E3_R] = RBI3(w_x, w_y, w_z, w_xy, w_xz, w_yz, w_xyz)

% Deviation from Additivity
DA = w_xyz - (w_x * w_y * w_z);
if DA <= 0 % synergistic
    normFactorDA = abs(0 - (w_x * w_y * w_z));
elseif DA > 0 % antagonistic
    normFactorDA = abs(min([w_x, w_y, w_z]) - (w_x * w_y * w_z));
end
DA_R = DA / normFactorDA;


% Accounting for pairwise interactions
% E3 = DA_xyz - (w_x * DA_yz) - (w_y * DA_xz) - (w_z * DA_xy); % emergent interaction
E3 = w_xyz - (w_x * w_yz) - (w_y * w_xz) - (w_z * w_xy) + (2 * w_x * w_y * w_z); % in terms of relative fitness
if E3 <= 0
    normFactor = abs(0 - (w_x * w_yz) - (w_y * w_xz) - (w_z * w_xy) + (2 * w_x * w_y * w_z));
    E3_R = E3 / normFactor;

elseif E3 > 0
    
    rescaleType = 3;
    switch rescaleType
    case 0 % Extension of two-drug rescaling method
        rescaleFactor = min([w_x, w_y, w_z]);
    case 1 % Buffering relative to pairwise drug effects
        rescaleFactor = min([w_xy, w_xz, w_yz]);
    case 2 % Buffering relative to single and pairwise drug effects
        rescaleFactor = min([w_x, w_y, w_z, w_xy, w_xz, w_yz]);
    case 3 % Emergent three-way interaction
        rescaleFactor = min([(w_x * w_yz), (w_y * w_xz), (w_z * w_xy)]);
    end
    normFactor = abs(rescaleFactor - (w_y * w_xz) - (w_z * w_xy) + (2 * w_x * w_y * w_z));

    % Case of antagonistic buffering
    if w_xyz <= rescaleFactor
        E3_R = E3 ./ normFactor;
                

    % Case of antagonistic suppression
    elseif w_xyz > rescaleFactor
        E3_R = 1 + (w_xyz - rescaleFactor) ./ abs(1 - rescaleFactor);
    end
    
    
    
    % not sure which the correct rescaling factor is, second one is presumably used by Beppler et al.
    % Tekin et al. 2016 describes the alternative rescaling methods
end


end % end of function