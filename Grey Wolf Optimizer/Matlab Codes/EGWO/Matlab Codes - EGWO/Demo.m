%% Demostrate how to invoke the EGWO function

% Please refer to the following paper on the details of the Enhanced Grey Wolf Optimizer (EGWO)  
% ---------------------------------------------------------------------------------------------
% Kaiping Luo,
% Enhanced grey wolf optimizer with a model for dynamically estimating the location of the prey
% Applied Soft Computing, 2019. 77: 225 - 235.
% Linkage: https://doi.org/10.1016/j.asoc.2019.01.025
% ---------------------------------------------------------------------------------------------
% Written by Kaiping Luo in Matlab R2018b.
% Copyright@2019: Beihang University.

%% Syntax
% fa = EGWO(alg,prob);
% [fa,xa,con] = EGWO(alg,prob);

% imput:
% alg: a struct for setting the parameters of the EGWO algorithm.
% alg.wps: a scalar for setting the wolf pack size.
% alg.max_nfe: a scalar for setting the maximal number of function evaluation.

% prob: a struct for setting the parameters of the problem to be solved.
% prob.dim: a scalar for setting the dimensional number of the problem to be solved.
% prob.ub: a horizontal vector for setting the upper bound of the problem to be solved.
% prob.lb: a horizontal vector for setting the lower bound of the problem to be solved.
% prob.fobj: a function handle for setting the objective function to be minized.

% output:
% fa: a scalar for recording the score of the position of the alpha wolf, namely, the minimal function value found
% xa: a horiaontal vector for recording the position of the alpha wolf, namely, the optimal solution found 
% con: a vertical vector for recording the best objective function value found at each iteration. 

%% 
clear 
clc

%% Example 1: 
disp('Demo 1: Solve the Sphere function using the EGWO algorithm.')

prob.dim = 30;
prob.fobj = @(x) sum(x.^2);
prob.lb = -100*ones(1,prob.dim);
prob.ub = 100*ones(1,prob.dim);  

alg.wps = 30;
alg.max_nfe = 10000*prob.dim;

[fa,xa,con] = EGWO(alg,prob);
sprintf('The minimal function value found: %e',fa)
figure
plot(con)
title('Convergence: Sphere')
xlabel('Number of function evaluation')
ylabel('Value of objective function')


%% Example 2:  
disp('Demo 2: Solve the benchmark function in the open CEC2017 test suite.')

seq = 2; % the sequence of the CEC2017 benchmark function, 1-30.
prob.dim = 30; % dim = 2,10,30,50,100. 
prob.fobj = @(x) cec17_func(x,seq);
prob.lb = -100*ones(1,prob.dim);
prob.ub = 100*ones(1,prob.dim);

alg.wps = 30;
alg.max_nfe = 10000*prob.dim;

[fa,xa,con] = EGWO(alg,prob);
sprintf('The minimal function value found: %e',fa)
figure
plot(con)
title('Convergence')
xlabel('Number of function evaluation')
ylabel('Value of objective function')
