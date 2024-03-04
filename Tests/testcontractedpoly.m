%% Test contractedpoly():
p = [1 2 3 4 0 6];
assert(all(contractedpoly(p,2,2) == [1 0 0; 2 1 0; 3 0 1; 4 2 0; 6 0 2], 'all'));

%% Test contractedpoly():
p = [1 2 3 1e-6 0 0];
assert(all(contractedpoly(p,2,2,1e-5) == [1 0 0; 2 1 0; 3 0 1], 'all'));

%% Test contractedpoly():
K = monomialsmatrix(3,3);
coef = 1:length(K);
A = [coef' K];
B = circshift(A,3);
assert(all(contractedpoly(expandedpoly(B),3,3) == A, 'all'))