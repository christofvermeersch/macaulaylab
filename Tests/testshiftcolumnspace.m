%% Test shiftnullspace() with random shift function/column-wise:
% Define the tolerance:
tol = 1e-10;
% Create a Macaulay matrix:
P = cell(2,1);
P{1} = [1 0 1; -1 2 0];
P{2} = [1 0 1; -2 1 0];
M = macaulay(systemstruct(P),2);

% Define the random shift polynomial:
G = [randn(2,1) eye(2)];
D = shiftcolumnspace(M,G,2,[1 2],[]);
assert(all(D{2} - [0; 2]< tol, 'all'))
assert(all(D{3} - [0; 4]< tol, 'all'))
