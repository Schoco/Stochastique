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

%% Parameters, play with them

N_tres = Np/3;

%% *** SEQUENTIAL MONTE CARLO METHOD ***

nx = 6;   % dimension of our particles
X = cell(Np,param.itmax +1);   % particles positions will be stored in X
XPoint = cell(Np,param.itmax +1);	% particles velocities will be stored in XPoint
XMean = cell(1, param.itmax +1);  % to store the mean of the predictions on position
XPointMean = cell(1, param.itmax +1);  % to store the mean of the predictions on velocity

% ** Step 1 : Initial prediction
for i = 1:Np
    X{i,1} = rand([2 param.N]).*40-20;  % the positions of the birds in each particle are uniformly distributed in the field of view
    XPoint{i,1} = zeros(2, param.N); % the initial velocities of each bird in each particle are 0
end

% ** Start of loop:
for k = 1:param.itmax+1 %for each time step
    
    % ** Step 2: Compute Prediction at time step k from the previous time step k-1 for each particle
    %Recomputation of the model
    if k > 1 %We can not predict at k = 1
        for i = 1:Np %for each particle
            for b = 1:param.N %for each bird in the particle
                norm_x = norm(X{i,k-1}(:,b));
                norm_x_point = norm(XPoint{i,k-1}(:,b));
                % Fhome
                if(0 <= norm_x && norm_x <= param.rf)
                    fhome=-4*param.rp/(param.rf^2)*norm_x*(norm_x-param.rf);
                else
                    fhome=0;
                end
                Fhome=-X{i,k-1}(:,b)*fhome;
                
                % Fvel
                Fvel=XPoint{i,k-1}(:,b)*param.vp*(1-norm_x_point/param.v0);
                
                % Finter
                Finter=[0;0];
                for inter=1:param.N %for each other bird
                    if (inter ~= b) %only for two distinct bird
                        dist=norm(X{i,k-1}(:,b)-X{i,k-1}(:,inter));
                        if(dist<param.df)
                            finter=-4*param.dp/((param.d0-param.df)^2)*(dist-param.d0)*(dist-param.df);
                        else
                            finter=0;
                        end
                        Finter=(X{i,k-1}(:,inter)-X{i,k-1}(:,b))*finter+Finter;
                    end
                end
                
                % Noise
                Fnoise=randn([2 1])*param.sigmaN^2;
                
                %Compute next step
                XPoint{i,k}(:,b) = param.ts*(Fhome+Fvel+Finter+Fnoise)+XPoint{i,k-1}(:,b);
                X{i,k}(:,b) = param.ts^2/2*(Fhome+Fvel+Finter)+param.ts*XPoint{i,k-1}(:,b)+Fnoise*param.ts^3/3+X{i,k-1}(:,b);
            end
        end
    end
    %At this point we have computed the prediction at time step k for the Np particles
    
    % ** Step 3 : Compute the weights   
    diff = zeros(2,param.N,Np);
    state_dist = zeros(1,Np); 
    weights = zeros(1,Np);
    for i=1:Np
        diff(:,:,i) = X{i,k} - cell2mat(y(:,k)'); % The difference between estimation and observation, as a 2x3, matrix
        state_dist(i) = norm(sqrt(sum(diff(:,:,i).^2,1))); %computes the distance between two states, as the norm of the vector composed of the euclidean distances between each of the birds
        weights(i) = normpdf(state_dist(i),0,sigmaObs); %we compute the weights
    end
    weights = weights ./ sum(weights); %Normalize ; Note: not needed when using randsample
    
    % ** Step 4 : Resample our particles using the weights
    ind_sample = randsample(Np,Np,true,weights);
    for i=1:Np
        X{i,k} = X{ind_sample(i),k};
        
        %if we use the prediction, uncomment this and comment the lines under:
        XPoint{i,k} = XPoint{ind_sample(i),k};
        
        %if we use the velocity estimator, uncomment this and comment the line above:
%         if k > 1
%             XPoint{i,k} = (cell2mat(y(:,k)') - cell2mat(y(:,k-1)'))./param.ts;
%         else
%             XPoint{i,k} = XPoint{ind_sample(i),k};
%         end
    end
    
    % ** Step 5 : Save the mean prediction after resampling, but before regularization
    Xsum = zeros(2,param.N);
    XPointsum = zeros(2,param.N);
    for i = 1:Np %for each particle
        Xsum = Xsum + X{i,k};
        XPointsum = XPointsum + XPoint{i,k};
    end
    XMean{k} = Xsum ./ Np;
    XPointMean{k} = XPointsum ./ Np;
    
    % ** Step 6 : Regularize by adding noise
    % If Neff < n/3, we use the regularized SMC
    if 1/sum(weights.^2) < N_tres
        epsilon = reshape(mvnrnd(zeros(nx,1), eye(nx), Np), 2, param.N, Np); % epsilon as defined in doc
        hopt = (4/(nx+2))^(1/(nx+4))*Np^(-1/(nx+4)); % optimal bandwidth of Gaussian kernel for dimension = 6
        Gamma = 1; %No whitening
        for i = 1:Np
            X{i,k} = X{i,k} + hopt*Gamma*epsilon(:,:,i); % add small error
        end
    end
    
end  % for t

% We want xf to be a N * itermax+1 cell array
xf = cell(param.N, param.itmax+1);
for i=1:param.itmax+1
   for n=1:param.N
      xf{n,i} =  XMean{i}(:,n);
   end
end

if(disp)
    main = figure(1);
    hold on;
    axis([-20 20 -20 20])
    plot(0,0,'or', 'MarkerFaceColor', 'r');
    h = zeros(3,1);
    for k=1:param.itmax
            for n=1:param.N
                h(n) = plot(xf{n,k}(1),xf{n,k}(2),'*k');
            end
        pause(10/100)
        delete(h);
    end
    hold off
    close(main)
end

end