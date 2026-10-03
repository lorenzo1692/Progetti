function scf_jacket = scf(param, type_cable, var, E_cbl_LTS, E_cbl_HTS)
%SCALING_FACTOR Compute scaling factor depending on cable type and parameter.
%
%   scf = scaling_factor(param, type_cable, var, E_cbl_LTS, E_cbl_HTS)
%
%   param        : input parameter
%   type_cable   : cell array of cable types ('LTS' or 'HTS')
%   var          : column index selecting the cable
%   E_cbl_LTS    : energy value for LTS cables
%   E_cbl_HTS    : energy value for HTS cables

    persistent p_LTS p_HTS xxx

    % Initialize persistent data only once
    if isempty(p_LTS)
        xxx = [1.087,1.136,1.190,1.250,1.316,1.389,1.471,1.563,1.667,...
               1.786,1.923,2.083,2.273,2.500,2.778];

        yyy_LTS = [1.01,1.03,1.06,1.10,1.16,1.22,1.29,1.36,1.43,1.50,...
                   1.57,1.64,1.71,1.78,1.85];

        yyy_HTS = [1.01,1.02,1.02,1.04,1.06,1.07,1.11,1.13,1.16,1.18,...
                   1.22,1.27,1.29,1.29,1.31];

        % Pre-fit the polynomials (degree 5)
        p_LTS = polyfit(xxx, yyy_LTS, 5);
        p_HTS = polyfit(xxx, yyy_HTS, 5);
    end

    % Default output if outside interpolation domain
    if param <= 1 || param >= 2.8
        scf = 1.5;
        return
    end

    % Select cable type
    isLTS = strcmp(type_cable{1, var}, 'LTS');

    % Choose polynomial
    if isLTS
        % E_cbl = E_cbl_LTS;  % Restore this if you need E_cbl elsewhere
        p = p_LTS;
    else
        % E_cbl = E_cbl_HTS;
        p = p_HTS;
    end

    % Evaluate scaling factor
    scf_jacket = polyval(p, param);

end
