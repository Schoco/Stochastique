function [xout] = generate_bird_flocks(param, disp)
x=zeros(2,param.itmax+1,param.N);
x_point=zeros(2,param.itmax+1,param.N);
%x_sec=zeros(2,param.itmax+1,param.N);
Finter=[0;0];
%init
x(2,1,:)=linspace(-5,5,param.N);
x(1,1,:)=linspace(3,3,param.N);
%loop
for i=2:param.itmax
   for j=1:param.N
       %Fhome
       if(norm(x(:,i,j))<=param.rf && norm(x(:,i,j))>=0)
         fhome=-4*param.rp/(param.rf^2)*norm(x(:,i,j))*(norm(x(:,i,j))-param.rf);
       else
          fhome=0;
       end
       Fhome=-x(:,i,j)*fhome;
       %Fvel
       Fvel=x_point(:,i,j)*param.vp*(1-norm(x_point(:,i,j))/param.v0);
       %Finter
       
       for inter=1:param.N
           dist=norm(x(:,i,j)-norm(x(:,i,inter)));
           if(dist<param.df)
             finter=-4*param.dp/((param.d0-param.df)^2)*(dist-param.d0)*(dist-param.df);
           else
             finter=0;
           end
           Finter=(x(:,i,inter)-x(:,i,j))*finter+Finter;
           
       end
       size(Finter)
       noiserecord=randn([2 1]);
       size(noiserecord)
       size(x_point(:,i,j))
       x_point(:,i+1,j)=param.ts*(Fhome+Fvel+Finter)+noiserecord*param.sigmaN^2*param.ts;
       x(:,i+1,j)=param.ts^2/2*(Fhome+Fvel+Finter)+param.ts*x_point(:,i+1,j)+noiserecord*param.sigmaN^2*param.ts^3/3+param.ts*x_point(:,i,j)+x(:,i,j); 
   end
   if(disp)
        
        for i2=1:param.itmax
            for j2=1:param.N
                figure(1);
                plot(0,0,'or');
                plot(x(1,i2,j2),x(2,i2,j2),'*');
                axis([-20 20 -20 20])
                hold on;
                
            end
            hold off
            pause(0.1)
        end
        xout=x;
   end
end

