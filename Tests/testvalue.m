
            % Define polynomials/affine points:
            p = {[1 0 0 0; 2 1 0 0; 3 0 1 1; 4 0 3 0; 5 2 2 2], ...
                [1 0 0; 1 1 0; 1 0 1; 1 2 0; 1 1 1; 1 0 2] ...
                [1 0; 2 1; 1 2]};
            system = systemstruct(p);
            x = {[1 -2 3], [0.5 2], 2};
            
            % Define the expected results in the Vandermonde basis:
            expected1 = {133, 8.75, 9};

            % Define the expected results in the Chebyshev basis:
            expected2 = {476, 11, 12};            

            for k = 1:3

                assert(evalmon(system.coef{k},system.supp{k},x{k}) == expected1{k})
            end

            for k = 1:3
                assert(evalcheb(system.coef{k},system.supp{k},x{k}) == expected2{k})
            end
  
            % Define matrix equations/affine points:
            P = {{[1 2; 3 4; 3 4]; [2 1; 0 1; 1 3]; zeros(3,2); ...
                zeros(3,2); [3 4; 2 1; 0 1]; [1 2; 4 2; 2 1];}, ...
                {ones(5,5); ones(5,5); ones(5,5); ones(5,5)}};
            n = {2, 3};
            dmax = {2, 1};
            x = {[2 2], [1 -2 3]};
%             
            % Define the expected results in the Vandermonde basis:
            expected1 = {[21 28; 27 18; 13 18], 3*ones(5,5)};

            % Define the expected results in the Chebyshev basis:
            expected2 = {[24 34; 39 24; 19 21] ,3*ones(5,5)};            
%     end
tol = 10e-12;

            for k = 1:2
            mep = mepstruct(P{k},dmax{k},n{k});
                assert(norm(evalmon(mep.coef{1},mep.supp{1},x{k}) - expected1{k}) < tol)
            end

               for k = 1:2
            mep = mepstruct(P{k},dmax{k},n{k});
                assert(norm(evalcheb(mep.coef{1},mep.supp{1},x{k}) - expected2{k}) < tol)
            end
% 
%     methods (Test,ParameterCombination = 'sequential')
%         function evalmepvander(testCase,P,n,dmax,x,expected1)
%             result = evalmepvander(P,n,dmax,x);
%             testCase.verifyEqual(result,expected1);
%         end
% 
%         function evalmep(testCase,P,n,dmax,x,expected1)
%             result = evalmep(P,n,dmax,x,@evalmepvander);
%             testCase.verifyEqual(result,expected1);
%         end
% 
%         function evalmepcheb(testCase,P,n,dmax,x,expected2)
%             result = evalmepcheb(P,n,dmax,x);
%             testCase.verifyEqual(result,expected2);
%         end
% 
%         function evalmepwithchebbasis(testCase,P,n,dmax,x,expected2)
%             result = evalmep(P,n,dmax,x,@evalmepcheb);
%             testCase.verifyEqual(result,expected2);
%         end
%     end
% 
% end