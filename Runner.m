% -----
% LINMA 1731 - Project
% Authors: ALEXE Simon & SCHOVAERS Corentin
% Date: 18 May 2018
% -----

function [xout] = Runner(RunnerNum)
%Runner runs both our simulation and the video example of one of the
%parameters set
%	RunnerNum one of the three parameter set (1, 2 or 3)
if(RunnerNum==1)
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
    disp=1;
    figure(2);
    load('video_1.mat');
    movie(gcf,Fr);
    xout = generate_bird_flocks(param,disp);
    
elseif(RunnerNum==2) 
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
    figure(2);
    load('video_2.mat')
    movie(gcf,Fr)

elseif(RunnerNum==3)
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
    figure(2);
    load('video_3.mat')
    movie(gcf,Fr)
	
else
	'Parameter RunnerNum has an unknown value.'
end