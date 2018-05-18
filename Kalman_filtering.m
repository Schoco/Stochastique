% -----
% LINMA 1731 - Project
% Authors: ALEXE Simon & SCHOVAERS Corentin
% Date: 18 May 2018
% -----

function [xf] = Kalman_filtering(param,Np,Yobs,sigmaObs,disp)
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
xf=cell(param.N,param.itmax);
%Covariance matrix of the process
CovarMatrixProcess=[param.sigmaN^2*param.ts^3/3 0 param.sigmaN^2*param.ts^2/2 0; 0 param.sigmaN^2*param.ts^3/3 0 param.sigmaN^2*param.ts^2/2;...
    param.sigmaN^2*param.ts^2/2 0 param.sigmaN^2*param.ts 0;0 param.sigmaN^2*param.ts^2/2 0 param.sigmaN^2*param.ts];
%Covariance matrix of the observation
RCovObser=[sigmaObs^2 0;0 sigmaObs^2];
%First Predicted covariance
PredictCov=eye(4);
%Initialize the estimator
x=linspace(1,1,param.itmax+1);
y=linspace(1,1,param.itmax+1);
xp=linspace(0,0,param.itmax+1);
yp=linspace(0,0,param.itmax+1);
x(1)=3;y(1)=rand*10-5;
xp=0.1;yp=0.1;

%Observation operation here identity  but we only observe x and y and not
%xpoint ypoint
H=[1 0 0 0; 0 1 0 0];
for j=1:param.N
for i=1:param.itmax
    %Predict
    %% Model use for the estimator
        norm_x = sqrt(x(i)^2+y(i)^2);
        norm_x_point = sqrt(xp(i)^2+yp(i)^2);
       % Fhome
       if(0 <= norm_x && norm_x <= param.rf)
         fhome=-4*param.rp/(param.rf^2)*norm_x*(norm_x-param.rf);
       else
          fhome=0;
       end
       Fhome=-[x(i) y(i)]*fhome;
       
	   % Fvel
       Fvel=[x(i) y(i)]*param.vp*(1-norm_x_point/param.v0);
       
	  
       
	   % Noise
	   Fnoise=randn([2 1])*param.sigmaN^2;
       %noiserecord=zeros([2 1]);
      
       %Compute next step

       Xpoint = param.ts*(Fhome+Fvel+Fnoise')+[xp(i) yp(i)];
       xp(i+1)=Xpoint(1);  %Bcos Matlab doesn't what double assignement
       yp(i+1)=Xpoint(2);
       Xpaspoint = param.ts^2/2*(Fhome+Fvel)+param.ts*[xp(i) yp(i)]+Fnoise'*param.ts^3/3+[x(i) y(i)]; 
       x(i+1)=Xpaspoint(1);
       y(i+1)=Xpaspoint(2);
     
    %% Jacobian values
    J1=1+param.ts^2/2*(4*param.rp/param.rf^2*(3*x(i+1)^2+y(i+1)^2-param.rf*(2*x(i+1)^2+y(i+1)^2)/sqrt(x(i+1)^2+y(i+1)^2)));
    J2=param.ts^2/2*4*param.rp/param.rf^2*(2*x(i+1)*y(i+1)-param.rf*(x(i+1)*y(i+1))/sqrt(x(i+1)^2+y(i+1)^2));
    J3=param.ts^2/2*param.vp*(1-(2*xp(i+1)^2+yp(i+1)^2)/(sqrt(xp(i+1)^2+yp(i+1)^2)*param.v0))+param.ts;
    J4=param.ts^2/2*param.vp*(-xp(i+1)*yp(i+1)/(param.v0*sqrt(xp(i+1)^2+yp(i+1)^2)));
    J5=param.ts^2/2*4*param.rp/param.rf^2*(2*x(i+1)*y(i+1)-param.rf*(x(i+1)*y(i+1))/sqrt(x(i+1)^2+y(i+1)^2));
    J6=1+param.ts^2/2*(4*param.rp/param.rf^2*(3*y(i+1)^2+x(i+1)^2-param.rf*(2*y(i+1)^2+x(i+1)^2)/sqrt(x(i+1)^2+y(i+1)^2)));
    J7=param.ts^2/2*param.vp*(-xp(i+1)*yp(i+1)/(param.v0*sqrt(xp(i+1)^2+yp(i+1)^2)));
    J8=param.ts^2/2*param.vp*(1-(2*yp(i+1)^2+xp(i+1)^2)/(sqrt(xp(i+1)^2+yp(i+1)^2)*param.v0))+param.ts;
    J9=4*param.ts*param.rp/param.rf^2*(3*x(i+1)^2+y(i+1)^2-param.rf*(2*x(i+1)^2+y(i+1)^2)/(sqrt(x(i+1)^2+y(i+1)^2)));
    J10=param.ts*4*param.rp/param.rf^2*(2*x(i+1)*y(i+1)-param.rf*x(i+1)*y(i+1)/sqrt(x(i+1)^2+y(i+1)^2));
    J11=1+param.ts*param.vp*(1-(2*xp(i+1)^2+yp(i+1)^2)/(param.v0*sqrt(xp(i+1)^2+yp(i+1)^2)));
    J12=param.ts*param.vp*(xp(i+1)*yp(i+1)/(param.v0*sqrt(xp(i+1)^2+yp(i+1)^2)));
    J13=param.ts*4*param.rp/param.rf^2*(2*x(i+1)*y(i+1)-param.rf*x(i+1)*y(i+1)/sqrt(x(i+1)^2+y(i+1)^2));
    J14=4*param.ts*param.rp/param.rf^2*(3*y(i+1)^2+x(i+1)^2-param.rf*(2*y(i+1)^2+x(i+1)^2)/(sqrt(x(i+1)^2+y(i+1)^2)));
    J15=param.ts*param.vp*(xp(i+1)*yp(i+1)/(param.v0*sqrt(xp(i+1)^2+yp(i+1)^2)));
    J16=1+param.ts*param.vp*(1-(2*yp(i+1)^2+xp(i+1)^2)/(param.v0*sqrt(xp(i+1)^2+yp(i+1)^2)));

    Jacobian=[J1 J2 J3 J4; J5 J6 J7 J8; J9 J10 J11 J12; J13 J14 J15 J16];
    PredictCov=Jacobian*PredictCov*Jacobian'+CovarMatrixProcess;
    
    %% Update the values of the estimator and the covariance estimation
    yTilde=(Yobs{1,i})'-[x(i+1) y(i+1)];
    ResiduCov=H*PredictCov*H'+RCovObser;
    NearOptKal=PredictCov*H'*inv(ResiduCov);
    Nawak=NearOptKal*yTilde';
    x(i+1)=Nawak(1);
    y(i+1)=Nawak(2);
    xp(i+1)=Nawak(3);
    yp(i+1)=Nawak(4);
    PredictCov=(eye(4)-NearOptKal*H)*PredictCov;
    %Transform to cell
    xf{j,i+1}=[x(i+1) y(i+1)];
end
xf{1,1}=[x(1) y(1)];
end
if(disp) %Display the estimator (black) and the true val (blue))
    
    for i2=1:param.itmax
            for j2=1:param.N
                main=figure(1);
                h(j2) = plot(xf{j2,i2}(1),xf{j2,i2}(2),'*k');hold on
                plot(Yobs{j2,i2}(1),Yobs{j2,i2}(2),'*b');
                plot(0,0,'or', 'MarkerFaceColor', 'r');
                axis([-20 20 -20 20])
            end
            hold off
        pause(10/100)
        
    end
    hold off
    close(main)
end

