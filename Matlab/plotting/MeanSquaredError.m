function [out] = MeanSquaredError(yhat,ydata,n,varargin)
%MEANSQUAREDERROR Calculate the Mean Squared Error for Two Arrays
%   This function takes in two arrays and calculates the mean-squared error
%   for those arrays. This is particularly useful for determining how well
%   an estimator predicts a dataset and giving the resulting error in the
%   same units as the input quantities (i.e. if the two arrays a & b have
%   units of Newtons, the standard error will have units of Newtons). The
%   standard error of regression is only used if the user passes "n" as a
%   variable argument---otherwise, the traditional "standard error" is
%   used. Mean Squared Error can be directly compared between datasets.
%
%   Note: size(yhat)==size(ydata) or an error is thrown.
%
%   a: DOUBLE ARRAY; the first array to use for the standard error.
%   
%   b: DOUBLE ARRAY; the second array to use for the standard error.
%   
%   n: INTEGER; the number of parameters used in the estimator model.
%
%   nvar: INTEGER; the user may also pass this to specify the size of
%   the running variance window used inside the standard error calculation.

if ~isempty(varargin)
    for i = 1:length(varargin)
        switch i
            case 1
                n = varargin{i};
            case 2
                nvar = varargin{i};
            otherwise
                error('You attempted to pass too many parameters to Standard Error. Please verify the function call.');
        end
    end
else
    n = 0;
    nvar = 3;
end

if size(yhat) ~= size(ydata)
    error('The arrays yhat and ydata given to StandardError() are NOT the same size. Please verify them and try again.');
end

if n < 0
    error('You have passed a negative value for the number of parameters in the estimator model, yhat. Please verify this value and try again.');
end

if n > 0
    % Mean Squared Error (MSE)
    out = sum(((yhat-ydata).^2)./(movvar(ydata,nvar).^2))./(length(ydata)-length(n));
else
    % Standard Error (S)
    out = sum(((yhat-ydata).^2)./(length(ydata)-1))./sqrt(length(ydata));
end

end