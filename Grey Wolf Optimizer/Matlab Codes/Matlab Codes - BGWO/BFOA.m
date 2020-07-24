function  [bf,nse,bx,con] = BFOA(prob,alg)
% ---------------------------------------------------------------------------------------------
% Rezoug, A., M. Bader-El-Den and D. Boughaci,
% Guided genetic algorithm for the multidimensional knapsack problem. Memetic Computing, 2017: p. 1--14.

% alg.N = 40;
% alg.S = 3;
% alg.L = 3;
% alg.b = 30;
% ---------------------------------------------------------------------------------------------
% Written by Kaiping Luo in Matlab R2018b.
% Copyright@2019: Beihang University.

%% initialization
Pg = 0.5*ones(1,prob.n);
F = zeros(alg.N,1);
bf = 0;
nse = 0;
Y = zeros(alg.S,prob.n);
Fy = zeros(alg.S,1);
con = zeros(alg.Max_nse,1);
h = 0;
while nse < alg.Max_nse
    X = generation_flies(alg.N,Pg);
    for i = 1:alg.N
        for j = 1:alg.S
            Y(j,:) = X(i,:);
            for l = 1:alg.L
                k = ceil(rand*prob.n);
                Y(j,k) = 1-Y(j,k);
            end
            [Y(j,:),Fy(j)] = Repair(Y(j,:),prob);
            nse = nse+1;            
            if nse >= alg.Max_nse
                return
            end
        end
        [f, b] = max(Fy);
        if f>F(i)
            F(i) = f;
            X(i,:) = Y(b,:);
            if f>bf
                bf = f;
                bx = Y(b,:);
                if ~isempty(prob.opt) && bf>=prob.opt
                    return
                end
            end            
        end
        con(h+1:nse) = bf;
        h = nse;
    end
    
    i1 = ceil(rand*alg.N);
    i2 = ceil(rand*alg.N);
    delt = bx + 0.5*(X(i1,:)-X(i2,:));
    Pg = 1./(1+exp(-alg.b*(delt-0.5)));
end



function X = generation_flies(N,P)
m = length(P);
X = zeros(N,m);
for i = 1:N
    for j = 1:m
        if rand < P(j)
            X(i,j) = 1;
        end
    end
end

function [x,f]= Repair(x,prob)
g = prob.R*(x')-prob.B;
f = x*prob.P;
for j = prob.Seq(end:-1:1)
    if any(g>0)
        if x(j)>0
            x(j) = 0;
            g = g-prob.R(:,j);
            f = f-prob.P(j);
        end
    end
end
for j = prob.Seq
    if x(j)==0 && all(g+prob.R(:,j)<=0)
        x(j) = 1;
        g = g+prob.R(:,j);
        f = f+prob.P(j);
    end
end