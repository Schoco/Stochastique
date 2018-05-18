% -----
% LINMA 1731 - Project
% Authors: ALEXE Simon & SCHOVAERS Corentin
% Date: 18 May 2018
% -----

% These is the 3rd proposed set of parameters.
param.N=3;
param.itmax=150;
param.ts=0.1;
param.rf=10;
param.rp=1.5;
param.v0=2;
param.vp=0.2;
param.d0=2;
param.df=3;
param.dp=0.2;
param.sigmaN=10;
disp=1;
xout = generate_bird_flocks(param,disp);