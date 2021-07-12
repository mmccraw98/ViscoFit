function sse = SSE_Voigt_Map(data,params,elasticSetting,fluidSetting)
    %SSE_Voigt Calculate the SSE for the Voigt model
    %   Calculate the Sum of Squared Errors for the Generalized
    %   Voigt Model according to the Lee and Radok indentation
    %   configuration, given a set of input parameters (params).
    %   This function performs fitting for all force curves
    %   separately, by taking in an additional index (idx) compared
    %   to the standard SSE_Voigt function. This is intended
    %   solely for Force Map analysis, wherein each pixel
    %   (containing a single force curve) is treated for analysis.
    %   The output is the same as for the SSE_Voigt function - a
    %   Sum of Squared Errors for that particular pixel.

    % Calculate test indentation
    test_indentations = LR_Voigt(params,data{1},data{2},data{3},data{5},data{6},data{7},elasticSetting,fluidSetting);

    switch data{7}
        case "spherical"
            beta = 1.5;
        case "conical"
            beta = 2;
    end

    % calculate global residual
    sse_global = sum(((data{4}.^beta)-test_indentations).^2);
    sse = sum(sse_global);

    [tauInds,modulusInds] = getParamIndices(params);
    ub = zeros(size(params))+eps;
    lb = zeros(size(params));

    ub(modulusInds) = 1e2;
    lb(modulusInds) = 1e-12;

    tauCenters = data{8}.*(10.^( (1:length(params(3:2:end)))-1 ));
    ub(tauInds) = tauCenters*10;
    lb(tauInds) = tauCenters/10;

    if length(params) > 1
        if fluidSetting
            ub(2) = 1;
            lb(2) = 0;
        end
    end

    if any(ub-params < 0) || any(params-lb < 0)
        sse = Inf;
    end
end % End Voigt SSE Map Function