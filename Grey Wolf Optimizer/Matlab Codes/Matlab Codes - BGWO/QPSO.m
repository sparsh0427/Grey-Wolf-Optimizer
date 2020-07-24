function [bf,nse,bx,con] = QPSO(prob,alg)
% ---------------------------------------------------------------------------------------------
% Haddar, B., et al., A hybrid quantum particle swarm optimization for the Multidimensional Knapsack Problem.
% Engineering Applications of Artificial Intelligence, 2016. 55(""): p. 1 - 13.

% Fields of the 'alg' struct
% alg.s = 20;
% alg.alpha = 0.1;
% alg.beta = 0.9;
% alg.c1 = 0.4;
% alg.c2 = 0.2;
% alg.c3 = 0.4;

% ---------------------------------------------------------------------------------------------
% Written by Kaiping Luo in Matlab R2018b.
% Copyright@2019: Beihang University.


%% random initialization
X = round(rand(alg.s,prob.n));
Y = zeros(alg.s,prob.n);
F = zeros(alg.s,1);
bf = 0;
nse = alg.s;
con = zeros(alg.Max_nse,1);
for i = 1:alg.s
    [X(i,:), F(i)] = repair(X(i,:),prob);
    if F(i)>bf
        bf = F(i);
        bx = X(i,:);
        if ~isempty(prob.opt) && bf>=prob.opt
            return
        end
    end
    con(i) = bf;
end
px = X;
pf = F;
h = nse;
while nse<alg.Max_nse 
    yg = alg.alpha*bx+alg.beta*(1-bx);
    for i = 1:alg.s
        yl = alg.alpha*px(i,:)+alg.beta*(1-px(i,:));
        Y(i,:) = alg.c1*Y(i,:)+alg.c2*yg+alg.c3*yl;
        for j = 1:prob.n
            if rand>Y(i,j)
                X(i,j) = 1;
            else
                X(i,j) = 0;
            end
        end
        [X(i,:),F(i)] = repair(X(i,:),prob);
        nse = nse+1;        
        if nse >= alg.Max_nse
            return
        end
        if F(i)>pf(i)
            pf(i) = F(i);
            px(i,:) = X(i,:);
        end
    end
    [f,b] = max(F);
    if f>bf
        bf = F(b);
        bx = X(b,:);
        [bx, bf, nse] = LocalSearch(bx,bf,prob,nse,alg.Max_nse);
    end
    con(h+1:nse) = bf;
    h = nse;
    if (~isempty(prob.opt) && bf>=prob.opt) || nse>=alg.Max_nse
        return
    end    
end


function [x, f] = repair(x,prob)
g = prob.R*x'-prob.B;
f = x*prob.P;
[x,f,g] = dropping(x,f,g,prob);
[x,f] = adding(x,f,g,prob);

function [x,f,nse] = LocalSearch(x,f,prob,nse,Max_nse)
y = x;
e = f;
while true
    for k = 1:prob.n
        z = y;
        h = e;
        if z(k)==1
            z(k) = 0;
            h = h-prob.P(k);
            g = prob.R*z'-prob.B;
            [z,h] = adding(z,h,g,prob);
        else
            z(k) = 1;
            h = h+prob.P(k);
            g = prob.R*z'-prob.B;
            if any(g>0)
                [z,h] = dropping(z,h,g,prob);
            end
        end
        if h>e
            e = h;
            y = z;
        end
        nse = nse+1;
        if nse >= Max_nse
            return
        end
    end
    if e>f
        f = e;
        x = y;
    else
        break
    end
end

function [x,f,g] = dropping(x,f,g,prob)
for j = prob.Seq(end:-1:1) %dropping
    if any(g>0)
        if x(j)>0
            x(j) = 0;
            g = g-prob.R(:,j);
            f = f-prob.P(j);
        end
    end
end

function [x,f,g] = adding(x,f,g,prob)
for j = prob.Seq % adding
    if x(j)==0 && all(g+prob.R(:,j)<=0)
        x(j) = 1;
        g = g+prob.R(:,j);
        f = f+prob.P(j);
    end
end

