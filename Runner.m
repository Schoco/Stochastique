%%Stochastique
RunnerNum=3;
disp=1;
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
    figure(2);
    load('video_1.mat');
    movie(gcf,Fr);
    generate_bird_flocks(param,disp);
    
elseif(RunnerNum==2) 
    Np=200;%Number of particles
    sigmaObs=1;%Corruption of observation
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
    figure(2);
    load('video_2.mat')
    movie(gcf,Fr)
    y=generate_bird_flocks(param,0);
    xf = Particle_filtering(param, Np, y, sigmaObs, disp);


else
    Np=200;%Number of particles
    sigmaObs=0.3;%Corruption of observation
    param.N=1;
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
    y=generate_bird_flocks(param,0);
    for i=1:param.itmax
        for j=1:param.N

        y{j,i}=y{j,i}(:)+randn([2 1])*sigmaObs^2;
        end
    end
    xf = Kalman_filtering(param,Np,y,sigmaObs,disp);
%     figure(2);
%     load('video_1.mat')
%     movie(gcf,Fr)
%     

end