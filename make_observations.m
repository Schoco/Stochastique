% -----
% LINMA 1731 - Project
% Authors: ALEXE Simon & SCHOVAERS Corentin
% Date: 18 May 2018
% -----

function [y] = make_observations(param, x, sigmaObs)
%Make_observations Coputes the corrupted observations of the position.
%   param is the same vector as defined in State model discretization
%   section
%   x is a cell of size param.N * (param.itmax + 1) containing vectors of
%   size 2*1 giving the real position of each bird, typically computed by
%   the function generate_bird_flocks
%   sigmaObs is the standard deviation of the noise to be added.

%init
y = cell(param.N, param.itmax + 1);

for t = 1:param.itmax+1
    for n = 1:param.N
        y{n,t}=x{n,t}+sigmaObs^2*randn([2 1]); %Add noise on real position to get the observed position
    end
end

end

