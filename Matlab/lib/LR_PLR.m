function [out] = LR_PLR(params,time,dt,indentation,tipSize,varargin)
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
if ~isempty(varargin)
    if length(varargin) > 1
        for i = 1:length(varargin)
            switch(i)
                case 1
                    nu = varargin{i};
                case 2
                    tipGeom = varargin{i};
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
        c = 4*tan(tipSize);
        beta = 2;
end

% Make our Dirac Delta array for the convolution portion
diracArray = zeros(size(time));
diracArray(time-dt<2*eps) = 1;

% Set t' equal to dt, per Efremov (2017)
t_prime = dt;

% Calculate Power Law Rheology Modulus
if length(params) > 1
    E = params(1).*( (1 + time./t_prime) .^ (-params(2)) );
else
    E = params(1).*diracArray./dt;
end

% Calculate the action integral quantity
convData = zeros(size(diracArray));
startInd = find(diracArray);
endInd = horzcat(find(diracArray)-1,numel(diracArray));
endInd(1) = [];
for i = 1:length(startInd)
    temp = convnfft(gradient(indentation(startInd(i):endInd(i)).^(beta)),E(startInd(i):endInd(i)),'full');
    convData(startInd(i):endInd(i)) = temp(1:(1+endInd(i)-startInd(i)));
end

% Trim the dataset to the region of interest, since the convolution gives
% an array that is length(indentation)+length(Q)+1, which is twice as long
% as our time array.
out = c.*convData.*dt;

end