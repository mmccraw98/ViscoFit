function [simParams] = getSimParams(condNumber)
%GETSIMPARAMS Get the parameters used in our viscoelastic simulation
%   Take in the simulation number, and output the Generalized Maxwell
%   viscoelastic parameters used for that test condition.
    simParams = [];
    switch condNumber
        case 1
            simParams = [1e4 1e4 1e-4];
        case 2
            simParams = [1e4 1e4 1e-4 1e4 1e-3];
        case 3
            simParams = [1e4 1e4 1e-4 1e4 1e-3 1e4 1e-2];
        case 4
            simParams = [1e4 1e4 1e-4 1e4 1e-3 1e4 1e-2 1e4 1e-1];
        case 5
            simParams = [5e5 1e4 1e-4];
        case 6
            simParams = [5e5 1e4 1e-4 1e4 1e-3];
        case 7
            simParams = [5e5 1e4 1e-4 1e4 1e-3 1e4 1e-2];
        case 8
            simParams = [5e5 1e4 1e-4 1e4 1e-3 1e4 1e-2 1e4 1e-1];
        case 9
            simParams = [5e7 1e4 1e-4];
        case 10
            simParams = [5e7 1e4 1e-4 1e4 1e-3];
        case 11
            simParams = [5e7 1e4 1e-4 1e4 1e-3 1e4 1e-2];
        case 12
            simParams = [5e7 1e4 1e-4 1e4 1e-3 1e4 1e-2 1e4 1e-1];
        case 13
            simParams = [1e4 5e5 1e-4];
        case 14
            simParams = [1e4 5e5 1e-4 5e5 1e-3];
        case 15
            simParams = [1e4 5e5 1e-4 5e5 1e-3 5e5 1e-2];
        case 16
            simParams = [1e4 5e5 1e-4 5e5 1e-3 5e5 1e-2 5e5 1e-1];
        case 17
            simParams = [5e5 5e5 1e-4];
        case 18
            simParams = [5e5 5e5 1e-4 5e5 1e-3];
        case 19
            simParams = [5e5 5e5 1e-4 5e5 1e-3 5e5 1e-2];
        case 20
            simParams = [5e5 5e5 1e-4 5e5 1e-3 5e5 1e-2 5e5 1e-1];
        case 21
            simParams = [5e7 5e5 1e-4];
        case 22
            simParams = [5e7 5e5 1e-4 5e5 1e-3];
        case 23
            simParams = [5e7 5e5 1e-4 5e5 1e-3 5e5 1e-2];
        case 24
            simParams = [5e7 5e5 1e-4 5e5 1e-3 5e5 1e-2 5e5 1e-1];
        case 25
            simParams = [1e4 5e7 1e-4];
        case 26
            simParams = [1e4 5e7 1e-4 5e7 1e-3];
        case 27
            simParams = [1e4 5e7 1e-4 5e7 1e-3 5e7 1e-2];
        case 28
            simParams = [1e4 5e7 1e-4 5e7 1e-3 5e7 1e-2 5e7 1e-1];
        case 29
            simParams = [5e5 5e7 1e-4];
        case 30
            simParams = [5e5 5e7 1e-4 5e7 1e-3];
        case 31
            simParams = [5e5 5e7 1e-4 5e7 1e-3 5e7 1e-2];
        case 32
            simParams = [5e5 5e7 1e-4 5e7 1e-3 5e7 1e-2 5e7 1e-1];
        case 33
            simParams = [5e7 5e7 1e-4];
        case 34
            simParams = [5e7 5e7 1e-4 5e7 1e-3];
        case 35
            simParams = [5e7 5e7 1e-4 5e7 1e-3 5e7 1e-2];
        case 36
            simParams = [5e7 5e7 1e-4 5e7 1e-3 5e7 1e-2 5e7 1e-1];
    end
end

