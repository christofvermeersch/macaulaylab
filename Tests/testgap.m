%% Test gap():
A = randn(1,1)*randn(1,5);
B = randn(3,2)*randn(2,5);
C = randn(10,3)*randn(3,5);
Z = [A; B; A; B; A ; A; C];
[dgap,ma] = gap(Z,3,3);
assert(dgap == 2)
assert(ma == 3);

%% Test gap():
A = randn(3,2)*randn(2,5);
B = randn(9,2)*randn(2,5);
C = randn(30,3)*randn(3,5);
Z = [A; B; A; B; A ; A; C];
[dgap,ma] = gap(Z,3,3,3,10e-10);
assert(dgap == 2)
assert(ma == 4);

%% Test gap():
A = randn(3,2)*randn(2,5);
B = randn(9,2)*randn(2,5);
C = randn(30,3)*randn(3,5);
Z = [A; B; A + 10e-6*randn(3,5); B; A ; A; C];
[dgap,ma] = gap(Z,3,3,3,10e-5);
assert(dgap == 2)
assert(ma == 4);

%% Test gap():
A = randn(3,2)*randn(2,5);
B = randn(9,2)*randn(2,5);
C = randn(30,3)*randn(3,5);
Z = [A; B; A + 10e-6*randn(3,5); B; A ; A; C];
[dgap,ma] = gap(Z,3,3,3,10e-10);
assert(dgap == 3)
assert(ma == 5);