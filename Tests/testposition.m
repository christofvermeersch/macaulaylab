%% Test position():
assert(position([0 0 0 0]) == 1)

%% Test position():
assert(position([1 1]) == 5)

%% Test position():
assert(all(position(monomialsmatrix(3,3)) == 1:20))

%% Test position():
assert(all(position([2 2 0; 0 0 0; 1 0 1]) == [24 1 7]))