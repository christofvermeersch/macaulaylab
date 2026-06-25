classdef testmepstruct < matlab.unittest.TestCase
    methods (Test)
        % Test the class constructor:
        function constructortest(testCase)
            A = randn(3,2);
            B = randn(3,2);
            C = randn(3,2);
            tensor = zeros(3,3,2);
            tensor(1,:,:) = A;
            tensor(2,:,:) = B;
            tensor(3,:,:) = C;
            mep = mepstruct({A,B,C},[0 0; 0 1; 1 0]);
            testCase.verifyClass(mep,'mepstruct');
            testCase.verifyInstanceOf(mep,'problemstruct');
            testCase.verifyEqual(mep.mat,{A,B,C});
            testCase.verifyEqual(mep.coef,{tensor});
            testCase.verifyEqual(mep.supp,{[0 0; 0 1; 1 0]});
            testCase.verifyEqual(mep.basis,NaN);
        end

        function constructorbasistest(testCase)
            A = [1 0; 0 1]; B = [0 1; 1 0];
            mep = mepstruct({A, B}, [0; 1], @chebyshev);
            testCase.verifyEqual(mep.basis, @chebyshev);
        end

        function propertiestest(testCase)
            A = ones(3,2); B = ones(3,2); C = ones(3,2);
            mep = mepstruct({A, B, C}, [0 0; 0 1; 1 0]);
            [s, dmax, m, di, nnze, k, l] = properties(mep);
            testCase.verifyEqual(s, 1);
            testCase.verifyEqual(dmax, 1);
            testCase.verifyEqual(m, 2);
            testCase.verifyEqual(di, {1});
            testCase.verifyEqual(k, 3);
            testCase.verifyEqual(l, 2);
        end
    end
end
