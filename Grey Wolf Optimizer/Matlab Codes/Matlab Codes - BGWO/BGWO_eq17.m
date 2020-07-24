function [fa,nse,xa] = BGWO_eq17(prob, alg)
% A binary grey wolf optimizer for the multidimensional knapsack problem
% using Equation 17
% ---------------------------------------------------------------------------------------------
% Kaiping Luo,Qiuhong Zhao.
% A binary grey wolf optimizer for the multidimensional knapsack problem.
% Applied Soft Computing, 2019. 83:
% Linkage: https://doi.org/10.1016/j.asoc.2019.105645
% ---------------------------------------------------------------------------------------------
% Written by Kaiping Luo in Matlab R2018b.
% Copyright@2019: Beihang University.


%% initialization
[WP,F,fa,fb,fc,xa,xb,xc,a,b,c] = initialization(prob,alg);
nse = alg.s;
x = zeros(1,prob.n);
%% genetic evolution
while nse < alg.Max_nse
    w = [fa fb fc]/(fa+fb+fc);
    u = exp(-100*nse/alg.Max_nse);
    xp = w*[xa;xb;xc]+u*randn(1,prob.n);     
    for i = 1:alg.s
        for j = 1:prob.n
            r = 2*(2*rand-1);
            y = xp(j) - r*abs(xp(j)-WP(i,j)); 
            r =  abs(y/sqrt(1+y^2)); % using equation 17
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
    for j = prob.Seq
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