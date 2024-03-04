%% Test nullrecrmacaulay():
tol = 1e-12;

M1 = macaulay(toy1,10);
M2 = macaulayupdate(M1,toy1,11);

Z = nullrecrmacaulay(null(full(M1)),full(M2(size(M1,1)+1:end,:)));
assert(norm(M2*Z) < tol);

%% Test nullsparsemacaulay():
tol = 1e-12;

M1 = macaulay(toy1,10);
M2 = macaulayupdate(M1,toy1,11);

Z = nullsparsemacaulay(null(full(M1)),toy1,11);
assert(norm(M2*Z) < tol);


%% Test nullrecrblockmacaulay():
tol = 1e-12;
P = cell(6,1);
for k = 1:6
    P{k} = randn(11,10);
end

M1 = macaulay(mepstruct(P,2,2),10);
M2 = macaulayupdate(M1,mepstruct(P,2,2),11);

Z = nullrecrmacaulay(null(full(M1)),full(M2(size(M1,1)+1:end,:)));
assert(norm(M2*Z) < tol);

%% Test nullsparseblockmacaulay():
tol = 1e-12;
P = cell(6,1);
for k = 1:6
    P{k} = randn(11,10);
end

M1 = macaulay(mepstruct(P,2,2),10);
M2 = macaulayupdate(M1,mepstruct(P,2,2),11);

Z = nullsparsemacaulay(null(full(M1)),mepstruct(P,2,2),11);
assert(norm(M2*Z) < tol);