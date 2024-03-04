A00 = randn(4,3);
A10 = randn(4,3);
A01 = randn(4,3);

p1 = [1 0 0; 2 1 0; 3 0 1];
p2 = [1 0 0; 1 2 2];

%% Test problemstruct():
problem = problemstruct({p1(:,1),p2(:,1)},{p1(:,2:3),p2(:,2:3)});
assert(all(problem.coef{1} == [1; 2; 3]))
assert(all(problem.coef{2} == [1; 1]))
assert(all(all(problem.supp{1} == [0 0; 1 0; 0 1])))
assert(all(all(problem.supp{2} == [0 0; 2 2])))
assert(problem.n == 2)
assert(problem.s == 2)
assert(problem.dmax == 4)
assert(problem.nnze{1} == 3)
assert(problem.nnze{2} == 2)
assert(problem.di{1} == 1)
assert(problem.di{2} == 4)

%% Test systemstruct():
system = systemstruct({p1,p2});
assert(all(system.coef{1} == [1; 2; 3]))
assert(all(system.coef{2} == [1; 1]))
assert(all(all(system.supp{1} == [0 0; 1 0; 0 1])))
assert(all(all(system.supp{2} == [0 0; 2 2])))
assert(system.n == 2)
assert(system.s == 2)
assert(system.dmax == 4)
assert(system.nnze{1} == 3)
assert(system.nnze{2} == 2)
assert(system.di{1} == 1)
assert(system.di{2} == 4)

%% Test mepstruct():
mep = mepstruct({A00, A10, A01},1,2);
tensor = mep.coef{1};
support = [0 0; 1 0; 0 1];
assert(all(all(squeeze(tensor(1,:,:)) == A00)))
assert(all(all(squeeze(tensor(2,:,:)) == A10)))
assert(all(all(squeeze(tensor(3,:,:)) == A01)))
assert(all(all(mep.supp{1} == support)))
assert(mep.n == 2)
assert(mep.s == 1)
assert(mep.dmax == 1)