function sse = SSE_PLR_Map(data,params,varargin)
    %SSE_PLR_Map Calculate the SSE for the PLR model
    %   Calculate the Sum of Squared Errors for the Power Law
    %   Rheology Model according to the Lee and Radok indentation
    %   configuration, given a set of input parameters (params).
    %   This function performs fitting for all force curves
    %   separately, by taking in an additional index (idx) compared
    %   to the standard SSE_PLR function above. This is intended
    %   solely for Force Map analysis, wherein each pixel
    %   (containing a single force curve) is treated for analysis.
    %   The output is the same as for the SSE_PLR function - a
    %   Sum of Squared Errors for that particular pixel.

    % Calculate test forces
    test_forces = LR_PLR(params,data{1},data{2},data{4},data{5},data{6},data{7});

    % calculate global residual
    sse_global = sum((data{3}-test_forces).^2);
    sse = sum(sse_global);

    % Power Law Rheology Roster:
    % [E_0 alpha]
    ub = [1e12;1];
    lb = [1e-2;0];

    if any(ub(1:length(params))-params < 0) || any(params-lb(1:length(params)) < 0)
        sse = Inf;
    end
end % End PLR SSE Map Function