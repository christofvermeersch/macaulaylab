classdef testnullity < matlab.unittest.TestCase
    methods (Test)
        function fullranktest(testCase)
            M = eye(3);
            testCase.verifyEqual(nullity(M), 0);
        end

        function rankdeficienttest(testCase)
            M = [1 2 3; 4 5 6; 7 8 9];
            testCase.verifyEqual(nullity(M), 1);
        end

        function zeronullitytest(testCase)
            M = magic(4);
            testCase.verifyEqual(nullity(M), size(M,2) - rank(M));
        end

        function tolerancetest(testCase)
            M = diag([1 1e-12 1]);
            testCase.verifyEqual(nullity(M, 1e-10), 1);
            testCase.verifyEqual(nullity(M, 1e-14), 0);
        end

        function rectangulartest(testCase)
            M = [1 0 0; 0 1 0];
            testCase.verifyEqual(nullity(M), 1);
        end
    end
end
