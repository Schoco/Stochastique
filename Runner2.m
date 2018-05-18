% -----
% LINMA 1731 - Project
% Authors: ALEXE Simon & SCHOVAERS Corentin
% Date: 18 May 2018
% -----

% These is the 2nd proposed set of parameters.
param.N=3;
param.itmax=150;
param.ts=0.1;
param.rf=15;
param.rp=3;
param.v0=2;
param.vp=1;
param.d0=3;
param.df=10;
param.dp=3;
param.sigmaN=2;
disp=1;
xout = generate_bird_flocks(param,disp);