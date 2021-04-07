function [sse] = SumOfSquaredResiduals(a,b)
%SUMOFSQUAREDRESIDUALS Calculate the Sum of Squared Residuals for Two 
%Arrays
%   This function takes in two arrays and calculates the sum of squared
%   residuals for those arrays. This is particularly useful for determining
%   how well an estimator predicts a dataset. Note: size(a)==size(b) or an 
%   error is thrown.
%
%   a: DOUBLE; the first array to use for the sum of squared residuals.
%   
%   b: DOUBLE; the second array to use for the sum of squared residuals.
sse = sum((a-b).^2,'all');
end