%___________________________________________________________________%
%  Grey Wold Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% Grey Wolf Optimizer
% function [Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
clear
clc
func_num = 1;
dim = 10;
lb = -100;
ub =  100;
SearchAgents_no = 10;
Max_iter = 50;
runs = 3;
fobj = str2func('cec14_func');

for iii=1:30
     func_num = iii;
    for jjj=1:runs
        iii,jjj
% SearchAgents_no = 30; % Number of search agents

% Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

% Max_iter=1000; % Maximum numbef of iterations

% Load details of the selected benchmark function
% [lb,ub,dim,fobj]=Get_Functions_details(Function_name);

% initialize alpha, beta, and delta_pos
Alpha_pos = zeros(1,dim);
Alpha_score = inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve = zeros(1,Max_iter);

iter = 1;% Loop counter
counter = 5;
flag =0;
% Main loop
while iter<Max_iter
    
    for i=1:size(Positions,1)  
        
       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub = Positions(i,:)>ub;
        Flag4lb = Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % Calculate objective function for each search agent
        fitness = feval(fobj, Positions(i,:)',func_num);
      %  fobj(Positions(i,:));
        
        % Update Alpha, Beta, and Delta
      if fitness < Alpha_score 
             Delta_score = Beta_score; % Update delta
            Delta_pos = Beta_pos ;
            
            Beta_score = Alpha_score; % Update beta
            Beta_pos = Alpha_pos;
            
            Alpha_score = fitness; % Update alpha
            Alpha_pos = Positions(i,:);
            flag = 1;
        
        elseif (fitness > Alpha_score && fitness < Beta_score)
             Delta_score = Beta_score; % Update delta
            Delta_pos = Beta_pos ;
            
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        
        
        elseif (fitness > Alpha_score && fitness > Beta_score && fitness < Delta_score)
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
            
         end
    end
    
    

      a = 2-iter*((2)/Max_iter); % a decreases linearly from 2 to 0

    
  
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
    end
    iter=iter+1;    
    Convergence_curve(iter) = Alpha_score;
end
       xbest(jjj,:) = Alpha_pos;
        fbest(iii,jjj)= Alpha_score;

    end
    f_mean(iii)=mean(fbest(iii,:));
end

% display(['The best solution obtained by GWO is : ', num2str(Alpha_pos)]);
% display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Alpha_score)]);
