%% Test macaulayupdate():
M1 = sort(macaulay(toy1,5),1);
M2 = sort(macaulayupdate(macaulay(toy1,4),toy1,5),1);
assert(sum(M1-M2,'all') == 0)

%% Test blockmacaulayupdate():
P = cell(6,1);
for k = 1:6
    P{k} = randn(11,10);
end

M1 = sort(macaulay(mepstruct(P,2,2),11),1);
M2 = sort(macaulayupdate(macaulay(mepstruct(P,2,2),10),mepstruct(P,2,2),11),1);
assert(sum(M1-M2,'all') == 0)