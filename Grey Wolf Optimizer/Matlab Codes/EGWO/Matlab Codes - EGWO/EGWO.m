function [fa,xa,con] = EGWO(alg,prob)
%% Enhanced Grey Wolf Optimizer (EGWO)
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

%% Initialization
dim = prob.dim;
fobj = prob.fobj;
if size(prob.lb,1)>size(prob.lb,2)
    lb = prob.lb';
else
    lb = prob.lb;
end
if size(prob.ub,1)>size(prob.ub,2)
    ub = prob.ub';
else
    ub = prob.ub;
end
wps = alg.wps;
max_nfe = alg.max_nfe;
clear prob alg % The struct data are readable but take more computational time.

fa = inf;
xa = zeros(1,dim);
fb = inf;
xb = zeros(1,dim);
fd = inf;
xd = zeros(1,dim);
X = zeros(wps,dim);
F = zeros(wps,1);
con = zeros(max_nfe,1);
nfe = 0;
for i = 1:wps
    X(i,:) = rand(1,dim).*(ub-lb)+lb;
    F(i) = fobj(X(i,:)');
    if F(i) < fa
        fa = F(i);
        xa = X(i,:);
    elseif F(i) > fa && F(i) < fb
        fb = F(i);
        xb = X(i,:);
    elseif F(i) > fb && F(i) < fd
        fd = F(i);
        xd = X(i,:);
    end
    nfe = nfe+1;
    con(nfe) = fa;
end
Max_iter = max_nfe/wps - 1;

%% Evaluation
for k = 1:Max_iter
    v = sort(rand(1,3),'descend');
    w = v/sum(v);
    u = exp(-100*k/Max_iter); % omega,the coefficient may be problem-dependence
    xp = w*[xa;xb;xd]+u*randn(1,dim);      
    for i = 1:wps
        for j = 1:dim
            r = 2*(2*rand-1); 
            y = xp(j) - r*abs(xp(j)-X(i,j));
            if y>ub(j)
                X(i,j) = X(i,j)+rand*(ub(j)-X(i,j));
            elseif y< lb(j)
                X(i,j) = X(i,j)+rand*(lb(j)-X(i,j));
            else
                X(i,j) = y;
            end
        end
        F(i) = fobj(X(i,:)');
        if F(i)<fa            
            fa = F(i);
            xa = X(i,:);
        elseif F(i)>fa && F(i)<fb
            fb = F(i);
            xb = X(i,:);
        elseif F(i)>fb && F(i)<fd
            fd = F(i);
            xd = X(i,:);
        end
    end
    con(k*wps+1:(k+1)*wps) = fa;
end