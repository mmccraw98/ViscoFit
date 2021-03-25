function [out] = LR_Maxwell(params,time,indentation,radius,varargin)
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
if ~isempty(varargin)
    for i = 1:length(varargin)
        switch(i)
            case 1
                nu = cell2mat(varargin{i});
        end
    end
end

% Calculate coefficient for the action integral
c = (8*sqrt(radius))/(3*(1-nu));

% Make our time matrix (for all the arms)
time_mat = row2mat(time,length(params(3:2:end)));

% Calculate the amount of relaxation that occurs for all arms in time. This
% quantity in time will be removed from the initial modulus response, Eg
Q_arms = sum(params(3:2:end)./params(4:2:end).*exp(-time_mat./params(4:2:end)),1);

% Add the initial modulus, Eg, such that it is relaxed to Ee as time
% stretches toward infinity
if params(2) ~= 0
    % Eg relaxes in time due to the steady-state viscosity stored inside
    % params(2)
    Q = sum(params(1:2:end)) - (params(1)./params(2).*exp(-time./params(2))) - Q_arms;
else
    Q = sum(params(1:2:end)) - Q_arms; % Relax Eg in time
end

% Calculate the action integral quantity
out = c.*convnfft(indentation.^(1.5),Q,'full').*(time(2)-time(1));

% Trim the dataset to the region of interest, since the convolution gives
% an array that is length(indentation)+length(Q)+1, which is twice as long
% as our time array.
out = out(1:length(time));

end