%% Test nbmonomials():
assert(nbmonomials(-1,3) == 0);

%% Test nbmonomials():
assert(nbmonomials(0,3) == 1);

%% Test nbmonomials():
assert(nbmonomials(1,3) == 4);

%% Test nbmonomials():
assert(nbmonomials(2,2) == 6);

%% Test nbmonomials():
assert(nbmonomials(3,3) == 20);

%% Test nbmonomials():
assert(nbmonomials(4,4) == 70);