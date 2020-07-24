function [bf,nse,bx,con] = HHS(prob,alg)
% ---------------------------------------------------------------------------------------------
% Zhang, B., et al., An effective hybrid harmony search-based algorithm for solving multidimensional knapsack problems.
% Applied Soft Computing, 2015. 29: p. 288 - 297.

% parameter settings:
% alg.hms = 40;
% alg.hmcr = 0.9; 0.74
% alg.parmax = 0.99;
% alg.parmin = 0.01;
% alg.gn = 20;
% alg.S = 2;
% alg.L = 3;
% alg.Max_nse = 100000;

% ---------------------------------------------------------------------------------------------
% Written by Kaiping Luo in Matlab R2018b.
% Copyright@2019: Beihang University.

%% random initialization
HM = round(rand(alg.hms,prob.n));
F = zeros(alg.hms,1);
bf = 0;
nse = 0;
con = zeros(alg.Max_nse,1);
for i = 1:alg.hms
    [HM(i,:),F(i)] = repair(HM(i,:),prob);
    nse = nse+1;
    if F(i) > bf
        bf = F(i);
        bx = HM(i,:);
        if ~isempty(prob.opt) && bf>=prob.opt
            return
        end
    end
    con(i) = bf;
end
t = 1;
tune = zeros(alg.gn,prob.n);
fv = zeros(alg.gn,1);
trail = zeros(alg.S,prob.n);
fu = zeros(alg.S,1);
%% iteration
while nse<alg.Max_nse
    par = alg.parmin+(alg.parmax-alg.parmin)*t/alg.Max_nse;
    for k = 1:alg.gn % gn =20
        a = ceil(rand*alg.hms);
        if rand<alg.hmcr %hmcr = 0.74
            for j = 1:prob.n
                tune(k,j) = HM(a,j);
                if rand<par
                    tune(k,j) = bx(j);
                end
            end
        else
            for j = 1:prob.n
                if rand<0.5
                    tune(k,j) = 1;
                else
                    tune(k,j) = 0;
                end
            end
        end
        [tune(k,:),fv(k)] = repair(tune(k,:),prob);
        nse = nse+1;
        con(nse) = bf;
        if nse >= alg.Max_nse
            return
        end
    end
    [fi,i] = max(fv);
    [wf,w] = min(F);
    if fi>wf
        F(w) = fi;
        HM(w,:) = tune(i,:);
        if fi>bf
            bf = fi;
            bx = tune(i,:);
            if ~isempty(prob.opt) && bf>=prob.opt
                return
            end
        end
    end
    
    for i = 1:alg.hms
        for k = 1:alg.S %S=2
            trail(k,:) = HM(i,:);
            for m = 1:alg.L % L=2
                b = ceil(rand*prob.n);
                trail(k,b) = 1-trail(k,b);
            end
            [trail(k,:), fu(k)]= repair(trail(k,:),prob);
            nse = nse+1;
            con(nse) = bf;
            if nse >= alg.Max_nse
                return
            end
        end
        for k = 1:alg.S
            if fu(k)>F(i)
                HM(i,:) = trail(k,:);
                F(i) = fu(k);
                if F(i)>bf
                    bf = F(i);
                    bx = HM(i,:);
                    if ~isempty(prob.opt) && bf>=prob.opt
                        return
                    end
                end
            end
        end
    end
    t = t+1;
end


function [x, f] = repair(x,prob)
%% repair the infeasible solution
g = prob.R*x'-prob.B;
f = x*prob.P;
for j = prob.Seq(end:-1:1)
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
for j = prob.Seq
    if x(j)==0 && all(g+prob.R(:,j)<=0)
        x(j) = 1;
        g = g+prob.R(:,j);
        f = f+prob.P(j);
    end
end