function sse = SSE_Maxwell_Map(data,params,elasticSetting,fluidSetting)
    %SSE_Maxwell Calculate the SSE for the Maxwell model
    %   Calculate the Sum of Squared Errors for the Generalized
    %   Maxwell Model according to the Lee and Radok indentation
    %   configuration, given a set of input parameters (params).
    %   This function performs fitting for all force curves
    %   separately, by taking in an additional index (idx) compared
    %   to the standard SSE_Maxwell function. This is intended
    %   solely for Force Map analysis, wherein each pixel
    %   (containing a single force curve) is treated for analysis.
    %   The output is the same as for the SSE_Maxwell function - a
    %   Sum of Squared Errors for that particular pixel.

    % Calculate test forces
    test_forces = LR_Maxwell(params,data{1},data{2},data{4},data{5},data{6},data{7},elasticSetting,fluidSetting);

    % calculate global residual
    sse_global = sum((data{3}-test_forces).^2);
    sse = sum(sse_global);

    [tauInds,modulusInds] = getParamIndices(params);
    ub = zeros(size(params))+eps;
    lb = zeros(size(params));

    ub(modulusInds) = 1e12;
    lb(modulusInds) = 1e-2;

    tauCenters = data{8}.*(10.^( (1:length(params(3:2:end)))-1 ));
    ub(tauInds) = tauCenters*10;
    lb(tauInds) = tauCenters/10;

    if length(params) > 1
        if fluidSetting
            ub(2) = max(tauCenters)*1e2;
            lb(2) = min(data{2});
        end
    end

    if any(ub-params < 0) || any(params-lb < 0)
        sse = Inf;
    end
end % End Maxwell SSE Map Function