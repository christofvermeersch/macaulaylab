%% Test gapstdmonomials();
c = [1 2 5 6];
[dgap,ma] = gapstdmonomials(c,3,3);
assert(ma == 4)
assert(dgap == 3)

%% Test gapstdmonomials();
c = [1 2 3 14];
[dgap,ma] = gapstdmonomials(c,3,3);
assert(ma == 3)
assert(dgap == 2)

%% Test gapstdmonomials();
c = [1 2 3 30];
[dgap,ma] = gapstdmonomials(c,3,3,2);
assert(ma == 3)
assert(dgap == 2)