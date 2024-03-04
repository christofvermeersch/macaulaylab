%% Test systeminfo():
P = cell(2,1);
P{1} = [1 2 2; 2 10 0];
P{2} = [2 0 0; 1 1 1];
system = systemstruct(P);
[s,n,di,dmax,nnze] = probleminfo(system);
assert(s == 2)
assert(n == 2)
assert(di{1} == 10)
assert(di{2} == 2)
assert(dmax == 10)
assert(nnze{1} == 2)
assert(nnze{2} == 2)

%% Test mepinfo():
P = cell(6,1);
P{1} = [1 2; 3 4; 3 4];
P{2} = [2 1; 0 1; 1 3];
P{3} = zeros(3,2);
P{4} = zeros(3,2);
P{5} = [3 4; 2 1; 0 1];
P{6} = [1 2; 4 2; 2 1];
[s,n,di,dmax,nnze] = probleminfo(mepstruct(P,2,2));
assert(s == 1)
assert(n == 2)
assert(di{1} == 2)
assert(dmax == 2)
assert(nnze{1} == 22)