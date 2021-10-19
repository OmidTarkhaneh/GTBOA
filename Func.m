
function z=Func(x,nvar,FuncNum)
%% 

  global NFE;
   if FuncNum==1
        z=0;
        NFE=NFE+1;
        for i=1:nvar
            % high conditioned elliptic function
%             z=z+(((10^6)^(i-1/nvar-1))*((x(i))^2));
             %Powell sum
             z=z+(abs(x(i)))^(i+1);
        end
  %%      
   elseif FuncNum==2
       NFE=NFE+1;
       z=0;
       for i=2:nvar
           z=z+((x(i))^2);
       end
       z=x(1)^2+10^6*z;
 %%      
   elseif FuncNum==3
       NFE=NFE+1;
       z=0;
       for i=2:nvar
           z=z+((x(i))^2);
       end
       z=(10^6)*(x(1)^2)+z;
  %%     
   elseif FuncNum==4
       NFE=NFE+1;
       z=0;
       for i=1:nvar-1
           z=z+(100*((x(i)^2-x(i+1))^2)+(x(i)-1)^2);
       end
 %%      
   elseif FuncNum==5
       NFE=NFE+1;
       z1=0;
       for i=1:nvar
           z1=z1+x(i)^2;
       end
       
       z2=0;
       for i=1:nvar
           z2=z2+cos(2*pi*x(i));
       end
       z=-20*exp(-0.2*sqrt((1/nvar)*z1))-exp((1/nvar)*z2)+20+exp(1);
 %%      
   elseif FuncNum==6
       NFE=NFE+1;
       kmax=20;
       a=0.5;
       b=3;
       z1=0;
       for i=1:nvar
           for k=0:kmax
               z1=z1+a^k*cos(2*pi*(b^k)*(x(i)+0.5));
           end
       end
       
       z2=0;
       for k=0:kmax
           z2=z2+a^k*cos(2*pi*(b^k)*0.5);
       end
       z=z1-nvar*z2;
  %%     
   elseif FuncNum==7
       NFE=NFE+1;
           z1=0;
           z2=1;
           for i=1:nvar
               z1=z1+((x(i)^2)/4000);
               z2=z2*cos(x(i)/sqrt(i));
           end
           z=z1-z2+1;
 %%          
   elseif FuncNum==8
       NFE=NFE+1;
       z=0;
       for i=1:nvar
           z=z+((x(i)^2)-(10*cos(2*pi*x(i)))+10);
           
       end
  %%     
   elseif FuncNum==9
       NFE=NFE+1;
       g=0;
       y=zeros(nvar,1);
       for i=1:nvar
           y(i)=x(i)+4.209687462275036*10^2;
           if abs(y(i))<=500
               g=g+y(i)*sin((abs(y(i)))^(1/2));
           elseif y(i)>500
               g=g+((500-mod(y(i),500))*sin(sqrt(abs(500-mod(y(i),500))))-(((y(i)-500)^2)/(10000*nvar)));
           elseif y(i)<-500
               g=g+((mod(abs(y(i)),500)-500)*sin(sqrt(abs(mod(abs(y(i)),500)-500)))-(((y(i)+500)^2)/(10000*nvar)));
           end
           z=418.9829*nvar-g;
       end
  %%     
   elseif FuncNum==10
       NFE=NFE+1;
       z2=1;
       for i=1:nvar
           z1=0;
           for j=1:32
               z1=z1+abs(2^j*x(i)-round(2^j*x(i)))/2^j;
           end
           z2=z2*(1+i*z1)^(10/(nvar^1.2));
       end
       z=10/nvar^2*z2-10/nvar^2;
 %%      
   elseif FuncNum==11
       NFE=NFE+1;
       z1=0;
       z2=0;
       for i=1:nvar
           z1=z1+(x(i)^2);
           z2=z2+x(i);
       end
       z=((abs(z1-nvar))^(1/4))+(0.5*z1+z2)/nvar+0.5;
  %%     
   elseif FuncNum==12
       NFE=NFE+1;
       z1=0;
       z2=0;
       for i=1:nvar
           z1=z1+(x(i)^2);
           z2=z2+x(i);
       end
       z=((abs((z1^2)-(z2^2)))^(1/2))+(0.5*z1+z2)/nvar+0.5;
 %%      
   elseif FuncNum==13
       NFE=NFE+1;
       z1=0;
       z2=1;
       y=zeros(nvar,1);
       y(nvar)=(100*(x(nvar)^2-x(1))^2 + (x(nvar)-1)^2);
       for i=1:nvar-1
           y(i)=(100*(x(i)^2-x(i+1))^2 + (x(i)-1)^2);
           z1=z1+y(i)^2/4000;
           z2=z2*cos(y(i)/sqrt(i));
       end
       z=z1-z2+1;
 %%      
   elseif FuncNum==14
       NFE=NFE+1;
       z1=0;
       for i=1:nvar-1
           z1=z1+0.5+(sin(sqrt((x(i)^2)+(x(i+1)^2)))^2-0.5)/(1+0.001*((x(i)^2)+(x(i+1)^2)))^2;
           
       end
       z=z1+0.5+(sin(sqrt((x(nvar)^2)+(x(1)^2)))^2-0.5)/(1+0.001*((x(nvar)^2)+(x(1)^2)))^2;
 %%      
   elseif FuncNum==15
       NFE=NFE+1;
       z1=0;
       k=10;
       for i=1:nvar
           z1=z1+cos(k*x(i))*exp(-x(i)^2/2);
       end
       z=1-((1/nvar)*z1);
%%       
   elseif FuncNum==16
       NFE=NFE+1;
       z=0;
       for i=1:nvar
           z=z+(x(i))^2;
       end
 %%      
   elseif FuncNum==17
       NFE=NFE+1;
        % Constants
       z1=0;
       z2=0;
       a=10;
       k=100;
       m=4;
       y=zeros(nvar,1);
       %calculting y
       for i=1:nvar
           y(i)=1+(x(i)+1)/4;
       end
       % calculating z1
       for i=1:nvar-1
           z1=z1+((y(i)-1)^2)*floor((1+10*(sin(pi*y(i+1))^2)));
       end
       % calculating z2
       for i=1:nvar
           if x(i)>a
               z2=z2+(k*(x(i)-a)^m);
           elseif x(i)<a && x(i)>-a
               z2=z2+0;
           elseif x(i)<-a
               z2=z2+(k*(-x(i)-a)^m);
           end
       end
       %final
       z=(pi/nvar)*(10*(sin(pi*y(1))^2)+z1+(y(nvar)-1)^2)+z2;
 %%      
   elseif FuncNum==18
       % Constants
       z1=0;
       z2=0;
       a=5;
       k=100;
       m=4;
              
       %calculating z1
       for i=1:nvar-1
           z1=z1+(x(i)-1)^2 *floor((1+sin(3*pi*x(i+1))^2));
       end
       %calculating z2
       for i=1:nvar
           if x(i)>a
               z2=z2+(k*(x(i)-a)^m);
           elseif x(i)<a && x(i)>-a
               z2=z2+0;
           elseif x(i)<-a
               z2=z2+(k*(-x(i)-a)^m);
           end
       end
       %final
       z=0.1*(sin(3*pi*x(1))^2+z1+(x(nvar)-1)^2*floor((1+sin(2*pi*x(nvar))^2)))+z2;
   %%
   elseif FuncNum==19
       z=0;
       for i=1:nvar
           z=z+i*x(i)^4+(rand-eps);
       end
%%       
   elseif FuncNum==20
       z=0;   
   for i=1:nvar
       for j=1:i
      z=z+(x(j)^2);
       end   
   end
 %%  
   elseif FuncNum==21
       max=x(1);
       for i=2:nvar
           if abs(x(i))>max
               max=abs(x(i));
           end
       end
       z=max;
%%       
   elseif FuncNum==22
       z1=0;
       z2=0;
       for i=1:nvar
           z1=z1+abs(x(i));
           z2=z2*abs(x(i));
       end
       z=z1+z2;
%%       
   elseif FuncNum==23
       z=0;
       for i=1:nvar
           z=z+((floor(x(i)+0.5))^2);
       end
%%       
   elseif FuncNum==24
       z=sum(abs(x.*sin(x)+0.1*x));
%%       
   elseif FuncNum==25
       z=0;
       for i=1:nvar
           z=z+x(i)^6*(2+sin(1/x(i)));
       end
%%       
   elseif FuncNum==26
       
        z=7*x(1)^2-6*sqrt(3)*x(1)*x(2)+13*x(2)^2;
       
%%       
   elseif FuncNum==27
       z=0;
       for i=2:nvar
           z=z+x(1)^2-x(1)*x(2)+x(2)^2;
       end
%%       
   elseif FuncNum==28
       z=0;
       for i=1:nvar
           z=z+((x(i)-1)^2)+(x(1)-x(i)^2)^2;
       end
%%       
   elseif FuncNum==29
       z=0;
       for i=1:nvar
           z=z+i*x(i)^2;
       end
%%       
   elseif FuncNum==30
       z=0;
       for i=1:nvar
           z=z+(floor(abs(x(i))));
       end
%%       
   elseif FuncNum==31
       z=0;
       for i=1:nvar
           z=z+418.9829-x(i)*sin(sqrt(abs(x(i))));
       end
%%       
   elseif FuncNum==32
%        z=0;
%        for i=1:nvar
%            z=z+(10^(i-1)*x(i)^2);
%        end

%2013 f169
       z=0;
       
       for i=1:nvar
           r=rand;
           z=z+r*(abs(x(i)))^i;
       end
%%       
   elseif FuncNum==33
%        z=0.5+(sin(sqrt((x(1)^2)+(x(2)^2)))^2-0.5)/(1+0.001*((x(1)^2)+(x(2)^2))).^2;
         z=0.5+(((sin(x(1)^2-x(2)^2))^2-0.5)/(1+0.001*(x(1)^2+x(2)^2)^2));
%%       
   elseif FuncNum==34
       z=x.*sign(x);
%%       
   elseif FuncNum==35
       z=0;
   for i=1:nvar
      z=z+floor(x(i));
   end
%%   
   elseif FuncNum==36
%%       
   elseif FuncNum==37
       z=abs((x(1)^2)+(x(2)^2)+(x(1)*x(2)))+abs(sin(x(1)))+abs(cos(x(2)));
%%       
   elseif FuncNum==38
       z=-200*exp(-0.2.*sqrt(x(1)^2+x(2)^2));
%%       
   elseif FuncNum==39
   elseif FuncNum==40
%%       
   elseif FuncNum==41
       z1=0;
       
       for i=2:nvar-2
           z1=z1+(x(i-1)+10*x(i))^2;
       end
       z=z1+5*(x(i+1)-x(i+2))^2+(x(i)-2*x(i+1))^4+10*(x(i-1)-x(i+2))^4;
    %%
   elseif FuncNum==42
       z=0;
       for i=1:nvar
           z=z+abs(x(i)^5-3*x(i)^4+4*x(i)^3+2*x(i)^2-10*x(i)-4);
       end
%%       
   elseif FuncNum==43
       z=0;
       for i=1:nvar
           z=z+(x(i)^2-i)^2;
       end
%%       
   elseif FuncNum==44
       z1=0;
       for i=1:nvar
           z1=z1+x(i)^2;
       end
       z=1-cos(2*pi*sqrt(z1))+0.1*(sqrt(z1));
%%       
   elseif FuncNum==45
       % 2013 f54
%        z1=0;
%        for i=1:nvar
%           z1=z1+ x(i)^2;
%        end
%        z=-exp(-0.5*z1);

     % 2013, F48
        z1=0;
        for i=2:nvar
        z1=z1+i*(2*x(i)^2-x(i-1))^2;
        end
        z=(x(1)-1)^2+z1;
       

%% 
   elseif FuncNum==46 
       z1=0;
       z2=0;
       for i=1:nvar
           z1=z1+cos(5*pi*x(i));
           z2=z2+x(i)^2;
       end
       z=-0.1*z1-z2;
      
   end
   
   
end