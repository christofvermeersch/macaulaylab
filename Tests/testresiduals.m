%% Test systemresiduals():
P = cell(2,1);
P{1} = [-1 2 0; 2 1 1; 1 0 2; 5 1 0; -3 0 1; -4 0 0];
P{2} = [1 2 0; 2 1 1; 1 0 2; -1 0 0];
system = systemstruct(P);

[maxresidual, allresiduals] = residuals(system,[1 0; 1 1; -1 2; 0 0]);
assert(maxresidual == 16)
assert(all(allresiduals == [0; 3; 16; norm([-4 -1])]))

%% Test mepresiduals():
tol = 1e-10;

P = cell(3,1);
P{1} = [1 1; 1 1; 1 1];
P{2} = [1 0; 1 0; 1 0];
P{3} = [0 2; 0 2; 0 2];
x = [2 1];
assert(residuals(mepstruct(P,2,1),x) < tol)

%% Test mepresiduals():
tol = 1e-10;

P = cell(3,1);
P{1} = [1 1; 1 1; 1 1];
P{2} = [1 0; 1 1; 1 0];
P{3} = [0 2; 0 2; 0 2];
x = [2 1];
assert(residuals(mepstruct(P,2,1),x) > tol)