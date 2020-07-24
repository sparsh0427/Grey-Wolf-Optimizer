function [fa,nse,xa,con] = BGWO(prob, alg)
% A binary grey wolf optimizer for the multidimensional knapsack problem
% ---------------------------------------------------------------------------------------------
% Kaiping Luo,Qiuhong Zhao.
% A binary grey wolf optimizer for the multidimensional knapsack problem.
% Applied Soft Computing, 2019. 83:
% Linkage: https://doi.org/10.1016/j.asoc.2019.105645
% ---------------------------------------------------------------------------------------------
% Written by Kaiping Luo in Matlab R2018b.
% Copyright@2019: Beihang University.

%% Syntax
% fa = BGWO(prob,alg);
% [fa,nse,xa,con] = BGWO(prob,alg);

% imput:
% alg: a struct for setting the parameters of the BGWO algorithm.
% alg.s: a scalar for setting the wolf pack size.
% alg.max_nse: a scalar for setting the maximal number of solution evaluation.

% prob: a struct for setting the parameters of the problem to be solved.
% prob.m: a scalar for setting the number of the knapsack constraints.
% prob.n: a scalar for setting the number of the optional items.
% prob.P: a vertical vector for setting the profit of items.
% prob.B: a vertical vector for setting the capacity constraints.
% prob.R: a matrix for setting the requirement of each item to each type of resource.
% prob.U: a horizental vector for recording the pseudo-utility of item.
% prob.Seq: a horizental vector for recording the descending order of the pseudo-utility of item

% prob.best: a scalar for recording the known optimal objection function value found by other algorithms.(Optional)
% prob.opt: a scalar for recording the optimal objection function value found by the exact solver. (Optinal)
% prob.x: a vertical vector for recording the exact solution (Optinal)
% prob.obj: a scalar for recording the optimal objection function value of the corresponding relaxation programming problem
        
% output:
% fa: a scalar, the score of the position of the alpha wolf, namely, the maximal function value found
% nse: a scalar, the number of solution evaluation.
% xa: a vector, the position of the alpha wolf, namely, the optimal solution found 
% con: a vector recording the best objective function value found during each solution evaluation. 


%% initialization
[WP,F,fa,fb,fc,xa,xb,xc,a,b,c] = initialization(prob,alg);
nse = alg.s;
x = zeros(1,prob.n);
con = zeros(alg.Max_nse,1);
con(1:nse) = fa;
%% evolution
while nse < alg.Max_nse
    w = [fa fb fc]/(fa+fb+fc);
    u = exp(-100*nse/alg.Max_nse); % omega,the coefficient may be problem-dependence
    xp = w*[xa;xb;xc]+u*randn(1,prob.n);     
    for i = 1:alg.s
        for j = 1:prob.n
            r = 2*(2*rand-1);
            y = xp(j) - r*abs(xp(j)-WP(i,j)); 
            r = abs(tanh(y));
            if rand < r
                x(j) = 1;
            else
                x(j) = 0;
            end
        end        
        [x, f] = Repair(x, prob);
        if f>F(i) 
            F(i) = f;
            WP(i,:) = x;
            
            if F(i)>fa
                fa = F(i);
                xa = WP(i,:);
                a = i;
                if ~isempty(prob.opt) && fa>=prob.opt
                    return
                end
                
            elseif F(i)<fa && F(i)>fb
                fb = F(i);
                xb = WP(i,:);
                b = i;
            elseif F(i)<fb && F(i)>fc
                fc = F(i);
                xc = WP(i,:);
                c = i;
            end
        
        elseif i~=a && i~=b && i~=c
            F(i) = f;
            WP(i,:) = x;
        end
        nse = nse+1;
        con(nse) = fa;
    end    
end


function [WP,F,fa,fb,fc,xa,xb,xc,a,b,c] = initialization(prob,alg)
WP = zeros(alg.s, prob.n);
F = zeros(alg.s, 1);
fa = 0;
fb = 0;
fc = 0;
xa = zeros(1,prob.n);
xb = xa;
xc = xa;
a = 1;
b = 1;
c = 1;
for i = 1:alg.s
    g = -prob.B;
    for j = prob.Seq %in descending order of the pseudo-utility
        if rand < 0.5 && max(g+prob.R(:,j)) <= 0
            WP(i,j) = 1;
            g = g+prob.R(:,j);
            F(i) = F(i)+prob.P(j);
        end
    end
    if F(i)>fa
        fa = F(i);
        xa = WP(i,:);
        a = i;
    elseif F(i)<fa && F(i)>fb
        fb = F(i);
        xb = WP(i,:);
        b = i;
    elseif F(i)<fb && F(i)>fc
        fc = F(i);
        xc = WP(i,:);
        c = i;
    end
end

function [x,f]= Repair(x,prob)
g = prob.R*(x')-prob.B;
f = x*prob.P;
for j = prob.Seq(end:-1:1) %in ascending order of the psuedo-utility
    if any(g>0)
        if x(j)>0
            x(j) = 0;
            g = g-prob.R(:,j);
            f = f-prob.P(j);
        end
    else
        break
    end
end
for j = prob.Seq % in descending order of the psuedo-utility
    if x(j)==0 && all(g+prob.R(:,j)<=0)
        x(j) = 1;
        g = g+prob.R(:,j);
        f = f+prob.P(j);
    end
end