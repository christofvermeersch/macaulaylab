%% Test nullrecrrow():
tol = 1e-12;
A = randn(10,5)*randn(5,10);
B = randn(10,3)*randn(3,10);
C = randn(10,2)*randn(2,10);

Z = nullrecrrow(null([A; B]), C);
assert(sum([A; B; C]*Z, 'all') < tol)

%% Test nullrecrrow():
tol = 1e-12;
A = randn(10,5)*randn(5,10);
B = randn(10,3)*randn(3,10);
C = randn(10,2)*randn(2,10);

Z = nullrecrrow(null([A; B; C]), A);
assert(sum([A; B; C; A]*Z, 'all') < tol)