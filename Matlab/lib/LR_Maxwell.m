function [out] = LR_Maxwell(params,time,dt,indentation,tipSize,varargin)
%LR_Maxwell Calculate Generalized Maxwell L&R Action Integral
%   This function generates the Generalized Maxwell Model response to an
%   applied spherical indentation (indentation) in time (time) with a probe 
%   of tip radius (radius). This model is an application of the Lee and
%   Radok (1960) indentation configuration, which is well-studied in the
%   viscoelastic parameter extraction literature; therein, a spherical
%   probe is assumed to indent a viscoelastic half-space that is previously
%   undisturbed.

% Check the varargin
nu = 0.5; % Default poisson's ratio of the sample to incompressible (0.5)
tipGeom = "spherical";
elasticSetting = 1;
fluidSetting = 0;
if ~isempty(varargin)
    if length(varargin) > 1
        for i = 1:length(varargin)
            switch(i)
                case 1
                    nu = varargin{i};
                case 2
                    tipGeom = varargin{i};
                case 3
                    elasticSetting = varargin{i};
                case 4
                    fluidSetting = varargin{i};
            end
        end
    else
        nu = varargin;
    end
end

% Calculate coefficient for the action integral
switch tipGeom
    case "spherical"
        c = (8*sqrt(tipSize))./(3*(1-nu));
        beta = 1.5;
    case "conical"
        c = (2.*tan(tipSize.*pi./180))./(pi.*(1-nu.^2));
        beta = 2;
end

% Make our Dirac Delta array for the elastic term
diracArray = zeros(size(time));
diracArray((time-dt)<2*eps) = 1;
    
if length(params) > 1
    % Make our time matrix (for all the arms)
    time_mat = row2mat(time,length(params(3:2:end)));

    % Calculate the amount of relaxation that occurs for all arms in time. This
    % quantity in time will be removed from the initial modulus response, Eg
    Q_arms = sum(params(3:2:end)./params(4:2:end).*exp(-time_mat./params(4:2:end)),1);

    % Add the initial modulus, Eg, such that it is relaxed to Ee as time
    % stretches toward infinity
    if fluidSetting
        % Eg relaxes in time due to the steady-state viscosity stored inside
        % params(2)
        Q = sum(params(1:2:end)).*(diracArray./dt) - (params(1)./params(2).*exp(-time./params(2))) - Q_arms;
    else
        Q = sum(params(1:2:end)).*(diracArray./dt) - Q_arms; % Relax Eg in time
    end
else
    Q = params(1).*(diracArray./dt);
end

% Calculate the action integral quantity
convData = zeros(size(diracArray));
startInd = find(diracArray);
endInd = horzcat(find(diracArray)-1,numel(diracArray));
endInd(1) = [];
for i = 1:length(startInd)
    temp = convnfft(indentation(startInd(i):endInd(i)).^(beta),Q(startInd(i):endInd(i)),'full');
    convData(startInd(i):endInd(i)) = temp(1:(1+endInd(i)-startInd(i)));
end

% Trim the dataset to the region of interest, since the convolution gives
% an array that is length(indentation)+length(Q)+1, which is twice as long
% as our time array.
out = c.*convData.*dt;

end