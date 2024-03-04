%% Test expandedpoly():
K = monomialsmatrix(3,3);
coef = 1:length(K);
A = [coef' K];
B = circshift(A,3);
assert(all(expandedpoly(B) == coef))