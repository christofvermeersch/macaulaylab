%% Test shiftnullspace() with random shift function/degree block-wise:
% Define the tolerance:
tol = 1e-10;
% Define a canonical basis of the null space:
Z = [1 1; 0 2; 0 4; 0 4; 0 8; 0 16];
% Define the random shift polynomial:
G = [randn(2,1) eye(2)];
D = shiftnullspace(Z,G,2,NaN);
assert(all(D{2} - [0; 2]< tol, 'all'))
assert(all(D{3} - [0; 4]< tol, 'all'))

%% Test shiftnullspace() with given shift function/degree block-wise:
% Define the tolerance:
tol = 1e-10;
% Define a canonical basis of the null space:
Z = [1 1; 0 2; 0 4; 0 4; 0 8; 0 16];
% Define the random shift polynomial:
G = [2 1 0; 3 0 1];
D = shiftnullspace(Z,G,2,NaN);
assert(all(D{1} - [0; 16]< tol, 'all'))

%% Test shiftnullspace() with random shift function/row-wise:
% Define the tolerance:
tol = 1e-10;
% Define a canonical basis of the null space:
Z = [1 1; 0 2; 0 4; 0 4; 0 8; 0 16];
% Define the random shift polynomial:
G = [randn(2,1) eye(2)];
D = shiftnullspace(Z,G,2,[1 2]);
assert(all(D{2} - [0; 2]< tol, 'all'))
assert(all(D{3} - [0; 4]< tol, 'all'))

%% Test shiftnullspace() with a user-specified polynomial basis:
% Define the tolerance:
tol = 1e-10;
% Define a canonical basis of the null space:
Z = [1 1; 0 2; 0 4; 0 4; 0 8; 0 16];
% Define the random shift polynomial:
G = [randn(2,1) eye(2)];
D = shiftnullspace(Z,G,2,NaN,1,@basismon);
assert(all(D{2} - [0; 2]< tol, 'all'))
assert(all(D{3} - [0; 4]< tol, 'all'))

