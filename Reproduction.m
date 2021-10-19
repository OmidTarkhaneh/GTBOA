function [ch1 ch2]=Reproduction(p1,p2,best,VarSize)

       x1=p1.Position;
       x2=p2.Position;
       y1=p1.Cost;
       y2=p2.Cost;
       best=best.Position;     
       
% Epsilon_max max drug effect-- epsilion=0 no drug efffect
% epsilion=1 high drug effect
% EC--- drug concentration
% C(t) refer to drug dose

       beta=unifrnd(0.1,0.5,VarSize);
       EC=abs(randn(size(x1))).^(beta);
       epsilon_max=rand();
       % Drug Efficiency
       epsilon=epsilon_max.*randn(size(x1))/(EC);
       
       % Drug efficiency positions
       ch1=(1-epsilon).*(best-x1);
       ch2=(1-epsilon).*(best-x2);
       alpha=rand(size(x1));
    
       m=alpha.*x1+(1-alpha).*(x2-ch1);
       n=alpha.*x2+(1-alpha).*(x1-ch2);
    
       ch1=m;
       ch2=n;
       
%        RF=(x1+x2)/rand;
%        ch1=x1+rand*(RF*randi(size(x1)));
%        ch2=x2+rand*(RF*randi(size(x2)));


end