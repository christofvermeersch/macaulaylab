classdef testtensorization < matlab.unittest.TestCase
    methods (Test)
        function singlematrixtest(testCase)
            A = [1 2; 3 4];
            tensor = tensorization({A});
            testCase.verifyEqual(size(tensor), [1 2 2]);
            testCase.verifyEqual(squeeze(tensor), A);
        end

        function twomatricestest(testCase)
            A0 = [1 2; 3 4];
            A1 = [5 6; 7 8];
            tensor = tensorization({A0, A1});
            testCase.verifyEqual(size(tensor), [2 2 2]);
            testCase.verifyEqual(squeeze(tensor(1,:,:)), A0);
            testCase.verifyEqual(squeeze(tensor(2,:,:)), A1);
        end

        function threematricestest(testCase)
            A0 = eye(3);
            A1 = 2*eye(3);
            A2 = 3*eye(3);
            tensor = tensorization({A0, A1, A2});
            testCase.verifyEqual(size(tensor), [3 3 3]);
            testCase.verifyEqual(squeeze(tensor(1,:,:)), A0);
            testCase.verifyEqual(squeeze(tensor(2,:,:)), A1);
            testCase.verifyEqual(squeeze(tensor(3,:,:)), A2);
        end

        function nonsquarematrixtest(testCase)
            A0 = [1 2 3; 4 5 6];
            A1 = [7 8 9; 0 1 2];
            tensor = tensorization({A0, A1});
            testCase.verifyEqual(size(tensor), [2 2 3]);
            testCase.verifyEqual(squeeze(tensor(1,:,:)), A0);
            testCase.verifyEqual(squeeze(tensor(2,:,:)), A1);
        end
    end
end
