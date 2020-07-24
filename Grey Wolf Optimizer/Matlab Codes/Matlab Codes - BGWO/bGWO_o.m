function [fa,nse,xa,con] = bGWO_o(prob, alg)
% ---------------------------------------------------------------------------------------------
% Emary, E., H.M. Zawbaa and A.E. Hassanien,
% Binary grey wolf optimization approaches for feature selection.
% Neurocomputing, 2016. 172: p. 371--381.

% ---------------------------------------------------------------------------------------------
% Written by Kaiping Luo in Matlab R2018b.
% Copyright@2019: Beihang University.


%% initialization
WP = round(rand(alg.s,prob.n));
F = zeros(alg.s,1);
fa = 0;
fb = 0;
fc = 0;
xa = zeros(1,prob.n);
xb = xa;
xc = xa;
nse = alg.s;
t = 1;
Max_ite = alg.Max_nse/alg.s;
con = zeros(alg.Max_nse,1);
%% genetic evolution
while true
    for i = 1:alg.s
        [WP(i,:), F(i)] = Repair(WP(i,:), prob);
        if F(i)>fa
            fa = F(i);
            xa = WP(i,:);
            if ~isempty(prob.opt) && fa>=prob.opt
                return
            end
        elseif F(i)<fa && F(i)>fb
            fb = F(i);
            xb = WP(i,:);
        elseif F(i)<fb && F(i)>fc
            fc = F(i);
            xc = WP(i,:);
        end
        con(nse-alg.s+i) = fa;
    end
    a = 2*(1-t/Max_ite);
    for i = 1:alg.s
        for j = 1:prob.n
            x1 = xa(j) - 2*a*(rand-1)*abs(2*rand*xa(j)-WP(i,j));             
            x2 = xb(j) - 2*a*(rand-1)*abs(2*rand*xb(j)-WP(i,j)); 
            x3 = xc(j) - 2*a*(rand-1)*abs(2*rand*xc(j)-WP(i,j)); 
            y = (x1+x2+x3)/3;
            r = 1/(1+exp(-10*(y-0.5)));
            if rand < r
                WP(i,j) = 1;
            else
                WP(i,j) = 0;
            end
        end   
        nse = nse+1; 
        if nse >= alg.Max_nse
            return
        end
    end 
    t = t+1;
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