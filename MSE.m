% -----
% LINMA 1731 - Project
% Authors: ALEXE Simon & SCHOVAERS Corentin
% Date: 18 May 2018
% -----

function [mse] = MSE(param, y, m)
%UNTITLED Summary of this function goes here
%   param is the same vector as defined in State model
%   y is the cell containing the correct positions
%   m is the cell containing the estimated positions

sumdist = 0;
for i = 1:param.N
    for k = 1:param.itmax
        sumdist = sumdist + norm(y{i,k}-m{i,k})^2;
    end
end
mse = param.N^-1 * param.itmax^-1 * sumdist;

end

