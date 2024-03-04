%% Test companion():
assert(all(companion({[3 1; 2 2; 1 3]}) == [0 0 0; 1 0 -3; 0 1 -2], 'all'))

%% Test companion():
tol = 1e-10;

P = cell(1,1);
P{1} = [-6 0; 1 4; -5 3; 5 2; 5 1];
assert(sum(sort(eig(companion(P))) - [-1; 1; 2; 3]) < tol)