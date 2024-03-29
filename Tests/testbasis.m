classdef testbasis < matlab.unittest.TestCase
    
    properties (TestParameter)
            % Define original and shift polynomial:
            A = {monomialsmatrix(3,3), monomialsmatrix(3,3), ...
                monomialsmatrix(3,3), monomialsmatrix(3,3), ...
                monomialsmatrix(3,3), monomialsmatrix(3,3), ...
                monomialsmatrix(3,3)};
            B = {[0 0 0], [1 0 0], [0 1 0], [0 0 1], [2 0 0], ...
                [0 1 1], [3 2 1]};
            C = {[1 1 1], [2 1 0], [1 1 0], [1 0 1], [2 0 0], ...
                [0 3 1], [1 2 3]};
            
            % Define the expected results in the Vandermonde basis:
            expected1 = {monomialsmatrix(3,3) + [0 0 0], ...
                monomialsmatrix(3,3) + [1 0 0], ...
                monomialsmatrix(3,3) + [0 1 0], ...
                monomialsmatrix(3,3) + [0 0 1], ...
                monomialsmatrix(3,3) + [2 0 0], ...
                monomialsmatrix(3,3) + [0 1 1], ...
                monomialsmatrix(3,3) + [3 2 1]};

            % Define the expected results in the Chebyshev basis:
            expected2 = {[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1], [3 1 0; 1 1 0; 3 1 0 ; 1 1 0; 3 1 0; 1 1 0; 3 1 0; 1 1 0], ...
                [1 2 0; 1 2 0; 1 0 0; 1 0 0; 1 2 0; 1 2 0; 1 0 0; 1 0 0], [1 0 2; 1 0 2; 1 0 2; 1 0 2; 1 0 0; 1 0 0; 1 0 0; 1 0 0], ...
                [4 0 0; 0 0 0; 4 0 0; 0 0 0; 4 0 0; 0 0 0; 4 0 0; 0 0 0], [0 4 2; 0 4 2; 0 2 2; 0 2 2; 0 4 0; 0 4 0; 0 2 0; 0 2 0], ...
                [4 4 4; 2 4 4;4 0 4; 2 0 4;4 4 2; 2 4 2;4 0 2; 2 0 2]};     
    end

    methods (Test,ParameterCombination = 'sequential')
        function testbasismon(testCase,A,B,expected1)
            [result, c] = basismon(A,B);
            testCase.verifyEqual(result,expected1);
            testCase.verifyEqual(c,1);
        end

        function testbasischeb(testCase,B,C,expected2)
            [result, c] = basischeb(B,C);
            testCase.verifyEqual(result,expected2);
            testCase.verifyEqual(c, ones(8,1)*0.125);
        end
    end
end