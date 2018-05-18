% -----
% LINMA 1731 - Project
% Authors: ALEXE Simon & SCHOVAERS Corentin
% Date: 18 May 2018
% -----

function [xf] = Particle_filtering(param, Np, y, sigmaObs, disp)
%Particle_filtering Monte Carlo particle filtering to estimate the position
%of the birds.
%   param is the same vector as defined in State model
%   Np is the number of particles in the filter (use 200 or 1000 as tyoucal
%   values)
%   y is a cell of size param.N x (param.itmax + 1) conaining vectors of
%   size 21 giving the observed position of each bird. These observations
%   results from the corruption of the true positions by additive white
%   Gaussian noise of standard deviation sigmaObs.
%   xf is a cell of size param.N x (param.itmax + 1) containing vectors of
%   size 2 x 1 giving the estimated (or filtered) position of each bird.

mu_w = 0; %Mean for the observed value noise
sqrt_sigmaObs = sqrt(sigmaObs); %Variance "" ""

out_noise_pdf = @(w) 1/sqrt((2*pi)^d_y*abs(det(sigmaObs))) * exp(-.5*(w-mu_w)'*inv(sigmaObs)*(w-mu_w));  % pdf of the output noise w_t

for i=1:param.itmax+1
    for j=1:param.N
        y{j,i}=y{j,i}+sigmaObs^2*randn([2 1]); %Add noise on real position to get the observed position
    end
end
% *** SEQUENTIAL MONTE CARLO METHOD ***

n = Np;   % sample set size. Sugg: 1e2
X = cell(n,param.itmax +1);   % particles will be stored in X
Xtilde = cell(n,param.itmax +1);  % to store the predictions
XPoint = cell(n,param.itmax +1);   % particles will be stored in X
XtildePoint = cell(n,param.itmax +1);  % to store the predictions


% ** Generate initial sample set {x_0^i,...,x_0^n}:

for i = 1:n
  X{i,1} = [3 (10*rand(1,1)-5)];  % Initialize the first prediction to a random value on the x=3 axis
end

% ** Start loop on time:

for t = 1:param.itmax+1
  
  % ** Prediction
    %Recomputation of the modell
  for i = 1:n
      norm_x = norm(X{i,t});
        norm_x_point = norm(XPoint{i,t});
       % Fhome
       if(0 <= norm_x && norm_x <= param.rf)
         fhome=-4*param.rp/(param.rf^2)*norm_x*(norm_x-param.rf);
       else
          fhome=0;
       end
       Fhome=-X{i,t}*fhome;
       
	   % Fvel
       Fvel=XPoint{i,t}*param.vp*(1-norm_x_point/param.v0);
       
	   % Finter
       Finter=[0;0];
       for inter=1:param.N
           if (inter ~= j)
               dist=norm(X{i,t}-X{inter,t});
               if(dist<param.df)
                 finter=-4*param.dp/((param.d0-param.df)^2)*(dist-param.d0)*(dist-param.df);
               else
                 finter=0;
               end
               Finter=(X{inter,t}-X{i,t})*finter+Finter;
           end
       end
       
	   % Noise
	   Fnoise=randn([2 1])*param.sigmaN^2;
       
       %Compute next step
       
     XtildePoint{i,t+1}=param.ts*(Fhome+Fvel+Finter+Fnoise)+XPoint{i,t};
    Xtilde{i,t+1} = param.ts^2/2*(Fhome+Fvel+Finter)+param.ts*XPoint{i,t}+Fnoise*param.ts^3/3+param.ts*XPoint{i,t}+X{i,t}; 
  end
  
  
  % ** Correction
  
  %y = y_true(:,t+1);  % y is the true output at time t+1
  
  weights = zeros(1,n);
  for i=1:n
    weights(i) = out_noise_pdf(y{}-(Xtilde{i,t+1}));  % HIDDEN
  end

  % Resample the particles according to the weights:
      % We are using Matlab
      ind_sample = randsample(n,n,true,weights);

  for i=1:n
    XPoint{i,t+1} = XtildePoint{ind_sample(i),t+1};
    X{i,t+1} = Xtilde{ind_sample(i),t+1};
  end
xf=X;
end  % for t

if(disp)
    main = figure(1);
    hold on;
    axis([-20 20 -20 20])
    plot(0,0,'or', 'MarkerFaceColor', 'r');
    h = zeros(3,1);
    for i2=1:param.itmax
            for j2=1:param.N
                h(j2) = plot(y{j2,i2},y(2,i2,j2),'*k');
            end
        pause(10/100)
        delete(h);
    end
    hold off
    close(main)
end




end
