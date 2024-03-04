%% Test macaulay():
eqs = cell(2,1);
eqs{1} = [1 0 0; 2 1 0; 3 0 1; 4 0 2];
eqs{2} = [5 0 0; 6 1 0];
system = systemstruct(eqs);

M = [1 2 3 0 0 4 0 0 0 0; 0 1 0 2 3 0 0 0 4 0; 0 0 1 0 2 3 0 0 0 4; ...
    5 6 0 0 0 0 0 0 0 0; 0 5 0 6 0 0 0 0 0 0; 0 0 5 0 6 0 0 0 0 0; ...
    0 0 0 5 0 0 6 0 0 0; 0 0 0 0 5 0 0 6 0 0; 0 0 0 0 0 5 0 0 6 0];

assert(sum(M-macaulay(system,3), 'all') == 0)

%% Test blockmacaulay():
P = cell(6,1);
for k = 1:6
    P{k} = randn(11,10);
end
Z = zeros(11,10);

M = [P{1} P{2} P{3} P{4} P{5} P{6} Z Z Z Z; Z P{1} Z P{2} P{3} Z ...
    P{4} P{5} P{6} Z; Z Z P{1} Z P{2} P{3} Z P{4} P{5} P{6}];

assert(sum(M-macaulay(mepstruct(P,2,2),3), 'all') == 0)


% %% Test blockmacaulayfast():
% P = cell(6,1);
% for k = 1:6
%     P{k} = randn(11,10);
% end
% 
% assert(sum(blockmacaulayfast(mepstruct(P,2,2),10) - ...
%     blockmacaulay(mepstruct(P,2,2),10,'fast'), 'all') == 0)

% %% Test macaulayoriginal() and macaulayfast():
% P = cell(2,1);
% P{1} = [1 0 0; 2 1 0; 3 0 1; 4 0 2];
% P{2} = [5 0 0; 6 1 0];
% system = systemstruct(P);
% 
% assert(sum(macaulayoriginal(system,3)-macaulay(system,3,'original'), 'all') == 0)
% assert(sum(macaulayfast(system,3)-macaulay(system,3,'fast'), 'all') == 0)
% 
% %% Test macaulaynoinitial();
% P = cell(2,1);
% P{1} = [1 0 0; 2 1 0; 3 0 1; 4 0 2];
% P{2} = [5 0 0; 6 1 0];
% system = systemstruct(P);
% 
% M1 = macaulay(system,3,'noinitial');
% M2 = macaulaynoinitial(system,3);
% 
% M = [1 2 3 0 0 4 0 0 0 0; 0 1 0 2 3 0 0 0 4 0; 0 0 1 0 2 3 0 0 0 4; ...
%     5 6 0 0 0 0 0 0 0 0; 0 5 0 6 0 0 0 0 0 0; 0 0 5 0 6 0 0 0 0 0];
% 
% assert(sum(M-M1, 'all') == 0)
% assert(sum(M-M2, 'all') == 0)
% 
% %% Test macaulaycheb():
% T = cell(2,1); 
% T{1} = [-1.5 0 0; 1 1 0; -1.5 0 2]; 
% T{2} = [2 1 0; -6 0 1]; 
% T = systemstruct(T);
% 
% M1 = macaulayfast(T,3,@basischeb);
% M2 = macaulay(T,3,'chebyshev');
% 
% M = [-1.5 1 0 0 0 -1.5 0 0 0 0; 0.5 -1.5 0 0.5 0 0 0 0 -1.5 0; ...
%     0 0 -2.25 0 1 0 0 0 0 -0.75; 0 2 -6 0 0 0 0 0 0 0; ...
%     1 0 0 1 -6 0 0 0 0 0; -3 0 0 0 2 -3 0 0 0 0; ...
%     0 1 0 0 0 0 1 -6 0 0; 0 -3 1 0 0 0 0 1 -3 0; ...
%     0 0 -3 0 0 0 0 0 2 -3];
% 
% assert(sum(M-M1, 'all') == 0)
% assert(sum(M-M2, 'all') == 0)