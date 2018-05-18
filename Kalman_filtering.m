% -----
% LINMA 1731 - Project
% Authors: ALEXE Simon & SCHOVAERS Corentin
% Date: 18 May 2018
% -----

function [xf] = Kalman_filtering(param,Np,y,sigmaObs,disp)
%Kalman_filtering Kalman filtering to estimate the position of the birds. 
%   param is the same vector as defined in State model
%   Np is the number of particles in the filter (use 200 or 1000 as tyoucal
%   values)
%   y is a cell of size param.N x (param.itmax + 1) conaining vectors of
%   size 21 giving the observed position of each bird. These observations
%   results from the corruption of the true positions by additive white
%   Gaussian noise of standard deviation sigmaObs.
%   xf is a cell of size param.N x (param.itmax + 1) containing vectors of
%   size 2 x 1 giving the estimated (or filtered) position of each bird.


end

