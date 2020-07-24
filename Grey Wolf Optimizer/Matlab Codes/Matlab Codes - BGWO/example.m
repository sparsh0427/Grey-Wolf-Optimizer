function example
% Demostrate how the procedure of the binary grey wolf optimizer is invoked

% Please refer to the following paper on the details of the binary grey wolf optimizer for the multidimensional knapsack problem
% ---------------------------------------------------------------------------------------------
% Kaiping Luo,Qiuhong Zhao.
% A binary grey wolf optimizer for the multidimensional knapsack problem.
% Applied Soft Computing, 2019. 83:
% Linkage: https://doi.org/10.1016/j.asoc.2019.105645
% ---------------------------------------------------------------------------------------------

% Written by Kaiping Luo in Matlab R2018b.
% Copyright@2019: Beihang University.

%% set the algorithmic parameters
clc
alg.s = 20;
alg.Max_nse = 100000;


%% Chu's benchmark problems
load('Benchmarks_Chu.mat') % 270 instances
k = 1 % 1-270
disp(['solve the Chu''s', int2str(k),'-th benchmark problem'])
prob = Problem(k)
f = BGWO(prob, alg) 
f13 = BGWO_eq13(prob, alg)
f14 = BGWO_eq14(prob, alg)
f15 = BGWO_eq15(prob, alg)
f17 = BGWO_eq17(prob, alg)
f18 = BGWO_eq18(prob, alg)
    
%% Weiing benchmark problems
load('Benchmarks_Weing.mat') % 8 instances
k = 2 % 1-8
disp(['solve the Weing', int2str(k),'-th benchmark problem'])
prob = Problem(k)
f = BGWO(prob, alg)
f13 = BGWO_eq13(prob, alg)
f14 = BGWO_eq14(prob, alg)
f15 = BGWO_eq15(prob, alg)
f17 = BGWO_eq17(prob, alg)
f18 = BGWO_eq18(prob, alg)

%% GK benchmark problems
load('Benchmarks_GK.mat') % 11 instances
k = 3 % 1-11
disp(['solve the GK', int2str(k),'-th benchmark problem'])
prob = Problem(k)
f = BGWO(prob, alg)  
f13 = BGWO_eq13(prob, alg)
f14 = BGWO_eq14(prob, alg)
f15 = BGWO_eq15(prob, alg)
f17 = BGWO_eq17(prob, alg)
f18 = BGWO_eq18(prob, alg)
 