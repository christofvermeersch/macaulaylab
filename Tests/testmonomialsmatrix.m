%% Test monomialsmatrix():
assert(all(monomialsmatrix(0,4) == [0 0 0 0]))

%% Test monomialsmatrix():
assert(all(monomialsmatrix(2,2) == [0 0; 1 0; 0 1; 2 0; 1 1; 0 2], 'all'))

%% Test monomialsmatrix():
assert(all(monomialsmatrix(1,3) == [0 0 0; 1 0 0; 0 1 0; 0 0 1], 'all'))