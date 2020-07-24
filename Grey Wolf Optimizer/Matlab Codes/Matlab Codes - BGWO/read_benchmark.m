function Problem=read_benchmark(Name)
% Function: 
% If Name="Chu", read Chu's MKP benchmark data including 270 problems. 
% If Name="GK", read GK MKP benchmark data including 11 problems.
% If Name="Weing", read Weing MKP benchmark data including 8 problems.
% If Name='all', read all MKP benchmark data
% Output Arguments:
% A data sturcts: Problem
% The field of Problem: m, n, P, R, B,o

% Written by Kaiping Luo in Matlab R2018b.
% Copyright@2019: Beihang University.

switch Name
    case 'GK'
        Problem=ReadGKData;
    case 'Weing'
        Problem=ReadWeingData;
    case 'Chu'
        Problem=ReadChuData;
    otherwise  %'all'
        ReadGKData;
        ReadWeingData;
        ReadChuData;
end

%%  Read the GK MKP Benchmarks
function Problem=ReadGKData
% 11 problems
mcell = cell(11,1);
ncell = cell(11,1);
pcell = cell(11,1);
rcell = cell(11,1);
bcell = cell(11,1);
kcell = cell(11,1);
scell = cell(11,1);
ucell = cell(11,1);
xcell = cell(11,1);
ocell = cell(11,1);
fcell = cell(11,1);
floder = pwd; % get the path of the current folder
for k=1:11
    file = [floder, '\GK MKP Benchmarks\gk', int2str(k),'.dat'];  
    v = dlmread(file, ' ', [0, 1,0,3]);
    ncell{k} = v(1);
    mcell{k} = v(2);
    kcell{k} = v(3);
    data1 = dlmread(file, ' ', 1, 1)';
    data2 = data1(data1 ~= 0);
    if length(data2)  ~= v(1) + v(2)*v(1) + v(2)
        error('It can not successfully read the data!')
    end
    pcell{k} = data2(1:v(1));
    rcell{k} = reshape(data2(v(1) +1: v(1) + v(2)*v(1)), [v(1), v(2)])';
    bcell{k} = data2(v(1) + v(1)*v(2) + 1:end);
    [scell{k},ucell{k},fcell{k}] = pseudo_utility(pcell{k},rcell{k},bcell{k});
    [xcell{k},f,flag] = Precise_Solution(-pcell{k},ncell{k},rcell{k},bcell{k});
    if flag == 1
        ocell{k} = -f;
    else
        ocell{k} = [];
    end
    fprintf('The %d-th GK MKP has been completed!\n',k)
end
Problem = struct('m',mcell,'n',ncell,'P', pcell,'R',rcell,'B',bcell,'opt',ocell,'x',xcell,'best',kcell,'obj',fcell,'Seq',scell,'U',ucell);
save('Benchmarks_GK.mat', 'Problem')

%%  Read the Weing MKP Benchmarks
function Problem=ReadWeingData
% 8 problems
mcell = cell(8,1);
ncell = cell(8,1);
pcell = cell(8,1);
rcell = cell(8,1);
bcell = cell(8,1);
kcell = cell(8,1);
scell = cell(8,1);
ucell = cell(8,1);
xcell = cell(8,1);
ocell = cell(8,1);
fcell = cell(8,1);
floder = pwd; % get the path of the current folder
for k=1:8
    file = [floder, '\Weing MKP Benchmarks\weing', int2str(k),'.dat'];  
    v = dlmread(file, ' ', [0,0,0,2]);
    ncell{k} = v(1);
    mcell{k} = v(2);
    kcell{k} = v(3);
    pcell{k} = dlmread(file, ' ', [1,0,1,v(1)-1])';
    rcell{k} = dlmread(file, ' ', [2,0,v(2)+1,v(1)-1]);
    bcell{k} = dlmread(file, ' ', [v(2)+2,0,v(2)+2,v(2)-1])';
    [scell{k},ucell{k},fcell{k}] = pseudo_utility(pcell{k},rcell{k},bcell{k}); 
    [xcell{k},f,flag] = Precise_Solution(-pcell{k},ncell{k},rcell{k},bcell{k});
    if flag == 1
        ocell{k} = -f;
    else
        ocell{k} = [];
    end
    fprintf('The %d-th Weing MKP has been completed!\n',k)
end
Problem = struct('m',mcell,'n',ncell,'P', pcell,'R',rcell,'B',bcell,'opt',ocell,'x',xcell,'best',kcell,'obj',fcell,'Seq',scell,'U',ucell);
save('Benchmarks_Weing.mat', 'Problem')

%%  Read the Chu's MKP Benchmarks
function Problem=ReadChuData
% 270 problems
% The number of constraints: m=[5, 10,30] 
% The number of variables: n=[100, 250, 500]
% The coefficient to generate the capacity for each resource: alpha=[0.25, 0.50, 0.75]
% The sequence in each group: s=1,2,10.
% The payoff coefficient: Pj=sum((Rij)i)/m +500qj, qj - U(0,1)
% The coefficient matrix: (Rij)m*n - U(0,1000)
% The capacity limitation: Bi=alpha*(sum(Rij)j), alpha=[0.25, 0.50, 0.75]
m = [5,10,30];
n  = [100, 250, 500];
alpha = [0.25, 0.50, 0.75];
mcell = cell(270,1);
ncell = cell(270,1);
acell = cell(270,1);
lcell = cell(270,1);
pcell = cell(270,1);
rcell = cell(270,1);
bcell = cell(270,1);
scell = cell(270,1);
ucell = cell(270,1);
xcell = cell(270,1);
ocell = cell(270,1);
fcell = cell(270,1);
floder = pwd; % get the path of the current folder

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:10
                h = 90*(i-1)+30*(j-1)+10*(k-1)+l;
                mcell{h} = m(i);
                ncell{h} = n(j);
                acell{h} = alpha(k);
                lcell{h} = l;
                
                file = [floder, '\Chu''s MKP Benchmarks\OR', int2str(m(i)), 'x', int2str(n(j)),'\OR', int2str(m(i)), 'x', int2str(n(j)),'-', sprintf('%0.2f',alpha(k)), '_', int2str(l),'.dat'];
                data1 = dlmread(file, ' ', 1, 1)';
                data2 = data1(data1 ~= 0);
                if length(data2)  ~= n(j) + m(i)*n(j) + m(i)
                    error('It can not successfully read the data!')
                end
                pcell{h} = data2(1:n(j));
                rcell{h} = reshape(data2(n(j) +1: n(j) + m(i)*n(j)), [n(j), m(i)])';
                bcell{h} = data2(n(j) + m(i)*n(j) + 1:end);
                
                [scell{h},ucell{h},fcell{h}] = pseudo_utility(pcell{h},rcell{h},bcell{h}); 
                [xcell{h},f,flag] = Precise_Solution(-pcell{h},ncell{h},rcell{h},bcell{h});
                if flag == 1
                    ocell{h} = -f;
                else
                    ocell{h} = [];
                end
                fprintf('The %d-th Chus MKP has been completed!\n',h)
            end 
        end
    end
end
file = [floder, '\Chu''s MKP Benchmarks\Best Known Values.xlsx'];                
kcell = mat2cell(xlsread(file,1,'B2:B271'),ones(270,1));
Problem = struct('m',mcell,'n',ncell,'a',acell, 's',lcell,'P', pcell,'R',rcell,'B',bcell,'opt',ocell,'x',xcell,'best',kcell,'obj',fcell,'Seq',scell,'U',ucell);
save('Benchmarks_Chu.mat', 'Problem')

function [seq,u,obj,pu] = pseudo_utility(P,R,B)
[m,n] = size(R);
f = [B; ones(n,1)];
A = -[R', eye(n)];
b = -P;
lb = zeros(m+n, 1);
ub = inf(m+n, 1);
options = optimoptions('linprog','Display','off');
[sprice,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);
% [sprice,~,exitflag] = cplexlp(f,A,b,[],[],lb,ub);
if exitflag == 1
    w = sprice(1:m);
    u = P'./(w'*R);
    [pu,seq] = sort(u,'descend');
    obj = fval;
else
    error('The dual model can not get its optimal solution!')
end

function [x,f,flag] = Precise_Solution(c,n,A,b)
options = optimoptions('intlinprog','Display','off','MaxTime',300);
[x,f,flag] = intlinprog(c,1:n,A,b,[],[],zeros(n,1),ones(n,1),[],options);