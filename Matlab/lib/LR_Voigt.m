function [out] = LR_Voigt(params,time,force,radius,varargin)
%LR_Voigt Calculate Generalized Kelvin-Voigt L&R Action Integral
%   This function generates the Generalized K-V Model response to the
%   force (F) that is applied by a spherical indenter in time (time) 
%   with tip radius (radius). This model is an application of the Lee and
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
c = (3*(1-nu))/(8*sqrt(radius));

% Make our time matrix (for all the arms)
time_mat = row2mat(time,length(params(3:2:end)));

% Calculate the amount of relaxation that occurs for all arms in time. This
% quantity in time will be removed from the initial modulus response, Eg
U_arms = sum(params(3:2:end)./params(4:2:end).*exp(-time_mat./params(4:2:end)),1);

% Add the initial compliance, Jg, such that it is increased to Je as time
% stretches toward infinity
if params(2) ~= 0
    % Steady-state fluidity is added to the Retardance
    U = params(1).*[1 zeros(1,length(time)-1)] + U_arms + params(2);
else
    U = params(1).*[1 zeros(1,length(time)-1)] + U_arms;
end

% Calculate the action integral quantity
out = c.*convnfft(force,U,'full').*(time(2)-time(1));

% Trim the dataset to the region of interest, since the convolution gives
% an array that is length(force)+length(Q)+1, which is twice as long as our
% time array.
out = out(1:length(time));

end