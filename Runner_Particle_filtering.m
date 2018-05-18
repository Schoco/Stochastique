% -----
% LINMA 1731 - Project
% Authors: ALEXE Simon & SCHOVAERS Corentin
% Date: 18 May 2018
% -----

% These is the 1st proposed set of parameters.
param.N=3;
param.itmax=150;
param.ts=0.1;
param.rf=10;
param.rp=3;
param.v0=2;
param.vp=2;
param.d0=0.8;
param.df=6;
param.dp=0.5;
param.sigmaN=0.2;
disp=0;

sigmaObs=1;
Np = 200;
runs = 15;

mse = zeros(1,runs);
for i = 1:runs
    x_real = generate_bird_flocks(param,disp);
    x_obs = make_observations(param, x_real, sigmaObs);
    x_estimated = Particle_filtering(param, Np, x_obs, sigmaObs, disp);
    mse(i) = MSE(param, x_real, x_estimated);
    i/runs
end

mse_mean = mean(mse);