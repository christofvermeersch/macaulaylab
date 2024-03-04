%% Test bezout():
P = cell(1,1);
P{1} = [1 5];
mb = bezout(systemstruct(P));
assert(mb == 5);

%% Test bezout():
P = cell(2,1);
P{1} = [1 3 0];
P{2} = [1 1 2];
mb = bezout(systemstruct(P));
assert(mb == 9);

%% Test bezout():
P = cell(3,1);
P{1} = [1 1 1 1];
P{2} = [1 0 3 0];
P{3} = [1 0 1 2];
mb = bezout(systemstruct(P));
assert(mb == 27);