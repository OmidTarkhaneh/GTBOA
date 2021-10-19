%___________________________________________________________________%
%  Golden Tortoise Beetle optimizer(GTBO) source codes demo version 1.0           %
%                                                                   %
%  Developed in MATLAB R2016a                                       %
%                                                                   %
%  Author and programmer: Omid Tarkhaneh                            %
%                                                                   %
%         e-Mail: TarkhanehOmid@gmail.com                           %
%                                                                   %
%                                                                   %
%                                                                   %
%                                                                   %
%   Main Code:                                                      %
%                                                                   %
%   Please just ran GTBO.m, to run different functions              %
%   simply change the value of 'F'                                  %
%                                                                   %
%                                                                   %
%___________________________________________________________________%







clc;
clear;
close all;
global NFE;


%% Problem Definition

emp_result.best=[];
emp_result.nfe=[];

TotalResult=repmat(emp_result,1,1);

kk=4;
for cc=kk:kk
    
F='F1';

[lb,ub,dim,fobj] = Get_Functions_details(F);  % Get the functions details


 xmin=lb;          % Lower Bound of Decision Variables
 xmax=ub;          % Upper Bound of Decision Variables
 
 nVar=dim;            % Number of Decision Variables

 
VarSize=[1 nVar];   % Decision Variables Matrix Size

%% DE Parameters

MaxIter=100;      % Maximum Number of Iterations
bestcost1=zeros(MaxIter,1);

nfe=zeros(MaxIter,2);
nPop=100;        % Population Size


%% Initialization
NFE=0;                        % Initial nvalue of Number of function evaluations
MaxNfe=100000;                % Maximum No of NFE

empty_individual.Position=[];
empty_individual.Cost=[];

fitness=10^10*ones(nPop,1);
globalcost=inf;
bestsol.Cost=inf;

MatureRate=0.8;                   %Mature Rate
SurvivalRate=0.2;                 % Survival Rate

SP=round(nPop*SurvivalRate);      % Number of Mature Beetle
IP=round(nPop*MatureRate);        % Number of Survival Larva


lamb_min=0.1;                          % Lower Bound of Initial WaveLength
lamb_max=0.9;                          % Upper Bound of Initial WaveLength

pop=repmat(empty_individual,nPop,1);
BestSol.Cost=inf; 

% Initial population of Beetles
for i=1:nPop
     
       pop(i).Position=xmin+(xmax-xmin).*rand(VarSize);
             
        pop(i).Cost=fobj(pop(i).Position);

        if pop(i).Cost<BestSol.Cost
            BestSol.Cost=pop(i).Cost;
            BestSol.Position=pop(i).Position;
        end

 end

BestCostBeetle=zeros(MaxIter,1);

%% Main Loop
it=1;
while NFE<MaxNfe

    
    SP=round(nPop*SurvivalRate);       % Number of Survival Larva
    IP=round(nPop*MatureRate);         % Number of Mature Beetle
      
    popi=repmat(empty_individual,IP,1);

    for i=1:IP
        
          s=pop(i).Position;
           
          A=randperm(nPop);
        
          A(A==i)=[];
        
          a=A(1);
          b=A(2);
          c=A(3);
           
          % Changin color using reflectivity and wavelength
           beta=unifrnd(lamb_min,lamb_max,VarSize);
           n_a=randn(VarSize).*beta;
           
           d_a=randn(VarSize);
           nb=randn(VarSize);
           
           a1 = 0.5+0.2*tan(pi*(rand(1,1)-1/2));
              if a1>1 || a1<0
               a1=0.5;
              end
           
           
           d_b=a1;
           
           step=n_a.*(d_a).*(cos(beta))+nb.*d_b.*cos(beta);
           
           wavelength(i)=(2.*rand().*sqrt(a1+sin(2*pi*rand()).^2))/rand();
           x=min([wavelength]);
           y=max([wavelength]);
           xx(i)=(wavelength(i)-x)/(y-x);
            if i==1
                xx(i)=0.5;
            end
            
            step=step-(rand()).*xx(i);
            
           %Generated Mutated Solution
           Infection=step.*(pop(a).Position-BestSol.Position);
        

          % Mutation  % Equation No 4.
            popi(i).Position=pop(i).Position+Infection;
            popi(i).Position=min(max(popi(i).Position , xmin),xmax);
            popi(i).Cost=fobj(popi(i).Position);
              
    end
    
     
%% Crossover for Survived Beetles

   pops=repmat(empty_individual,round(SP),2);
   for w=1:SP
%       p1=BinaryTournamentSelection(pop);
%       p2=BinaryTournamentSelection(pop);
     
    % Random selection
         A=randperm(nPop);
        
          A(A==i)=[];
        
          p1=A(1);
          p2=A(2);
          c=A(3);
      
      p1=pop(p1);
      p2=pop(p2);
      
      % Crossover strategy Equation No 7.
      [pops(w,1).Position, pops(w,2).Position]=Reproduction(p1,p2,BestSol,VarSize);
       pops(w,1).Position=min(max(pops(w,1).Position , xmin),xmax);
       pops(w,2).Position=min(max(pops(w,2).Position , xmin),xmax);
       pops(w,1).Cost=fobj(pops(w,1).Position);
       pops(w,2).Cost=fobj(pops(w,2).Position);
   end
   pops=pops(:);
   %% Update pop (Sort + remove extras)
   
    pop=[pop
         popi
         pops];
    Costs=[pop.Cost];
    [Costs SortOrder]=sort(Costs);
    pop=pop(SortOrder);
    
     % Removing extra Individuals
    pop=pop(1:nPop);
    Costs=Costs(1:nPop);

    % Saving the Results
    BestSol=pop(1);
    BestCostBeetle(it)=Costs(1);
    nfe(it,1)=NFE;
    nfe(it,2)=Costs(1);
     disp(['iteration= ',num2str(it),' bestcost= ',num2str(BestCostBeetle(it)),'func eva:=' num2str(nfe(it))]);
it=it+1;
end
 
BestRepCR(1)=BestCostBeetle(end);
medianCR=median(BestRepCR);
meanCR=mean(BestRepCR);
BestOfAllCR=min(BestRepCR);
WorstCR=max(BestRepCR);
StdCR=std(BestRepCR);
disp([ 'The Best CR = ' num2str(BestOfAllCR)]);
disp([ 'Mean CR = ' num2str(meanCR)]);
disp([ 'Median CR = ' num2str(medianCR)]);
disp([ 'Worst CR = ' num2str(WorstCR)]);
disp([ 'Std CR = ' num2str(StdCR)]);
 
 fprintf('next run  %s .\n ' ,'*****************************');

TotalResult(cc).best=BestOfAllCR;
TotalResult(cc).nfe=nfe;
save (['GTBOAResults_' num2str(cc) ',mat'], 'TotalResult', 'nfe','BestOfAllCR');
end






