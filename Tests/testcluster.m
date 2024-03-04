%% Test cluster():
tol = 1e-3;
D = cell(1,1);
cluster1 = 3*ones(4,1) + 1e-4*randn(4,1);
cluster2 = 2*ones(3,1) + 1e-4*randn(3,1);
cluster3 = ones(1,1);
D{1} = [cluster1; cluster2; cluster3];

[~,centers,clusters] = cluster(D,tol);
assert(sum(abs(centers - [3; 2; 1])) < 3*tol)
assert(all(clusters{1} == [1; 2; 3; 4]))
assert(all(clusters{2} == [5; 6; 7]))
assert(all(clusters{3} == 8))

%% Test cluster():
tol = 1e-3;
D = cell(2,1);
c1 = 3*ones(1,1) + 1e-4*randn(1,1);
c2 = 2*ones(1,1) + 1e-4*randn(1,1);
c3 = ones(1,1);
D{1} = [c1; c2; c3; c1; c1; c2; c1; c2];
D{2} = randn(8,1);

[E,centers,clusters] = cluster(D,tol);
assert(length(centers) == 3)
cl = zeros(length(centers),1);
for k = 1:length(centers)
    cl(k) = length(clusters{k});
end
[~,idx] = sort(cl);
assert(all(clusters{idx(3)} == [1; 4; 5; 7]))
assert(all(clusters{idx(2)} == [2; 6; 8]))
assert(all(clusters{idx(1)} == 3))