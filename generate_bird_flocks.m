function [xout] = generate_bird_flocks(param, disp)
x=zeros(2,param.itmax+1,param.N);
x_point=zeros(2,param.itmax+1,param.N);
%x_sec=zeros(2,param.itmax+1,param.N);
%init
x(2,1,:)=linspace(-5,5,param.N);
x(1,1,:)=linspace(3,3,param.N);
%loop
for i=1:param.itmax
   for j=1:param.N
		norm_x = norm(x(:,i,j));
        norm_x_point = norm(x_point(:,i,j));
       % Fhome
       if(0 <= norm_x && norm_x <= param.rf)
         fhome=-4*param.rp/(param.rf^2)*norm_x*(norm_x-param.rf);
       else
          fhome=0;
       end
       Fhome=-x(:,i,j)*fhome;
       
	   % Fvel
       Fvel=x_point(:,i,j)*param.vp*(1-norm_x_point/param.v0);
       
	   % Finter
       Finter=[0;0];
       for inter=1:param.N
           if (inter ~= j)
               dist=norm(x(:,i,j)-x(:,i,inter));
               if(dist<param.df)
                 finter=-4*param.dp/((param.d0-param.df)^2)*(dist-param.d0)*(dist-param.df);
               else
                 finter=0;
               end
               Finter=(x(:,i,inter)-x(:,i,j))*finter+Finter;
           end
       end
       
	   % Noise
	   Fnoise=randn([2 1])*param.sigmaN^2;
       %noiserecord=zeros([2 1]);
       %size(noiserecord);
       %size(x_point(:,i,j));
       
       %Compute next step
       x_point(:,i+1,j)=param.ts*(Fhome+Fvel+Finter+Fnoise)+x_point(:,i,j);
       x(:,i+1,j)=param.ts^2/2*(Fhome+Fvel+Finter)+param.ts*x_point(:,i+1,j)+Fnoise*param.ts^3/3+param.ts*x_point(:,i,j)+x(:,i,j); 
   end
end

if(disp)
    main = figure(1);
    hold on;
    axis([-20 20 -20 20])
    plot(0,0,'or', 'MarkerFaceColor', 'r');
    h = zeros(3,1);
    for i2=1:param.itmax
            for j2=1:param.N
                h(j2) = plot(x(1,i2,j2),x(2,i2,j2),'*k');
            end
        pause(5/100)
        delete(h);
    end
    hold off
    close(main)
end
xout=x;
end