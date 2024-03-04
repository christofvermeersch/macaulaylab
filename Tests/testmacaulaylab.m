classdef testmacaulaylab < matlab.unittest.TestCase

    properties (TestParameter)
            % Define problems:
            systems = {dreesen1, dreesen2, dreesen3, dreesen4, dreesen5};
%             systemscheb = {cheb1, cheb2, cheb3}
            blocked = {true, false};
            recursive1 = {'recursive', 'iterative','sparse'};
            recursive2 = {'recursive', 'iterative'};
            l = {5};
            n = {1, 2};
            dmax = {1, 2};
    end

    methods (Test)
        function solvesystemwithnull(testCase,systems,blocked,recursive1)
            options.blocked = blocked;
            options.recursive = recursive1;
            options.algorithm = 'null';
            [~,output] = macaulaylab(systems,20,options);
            result = output.accuracy;
            testCase.verifyLessThan(result,10e-10);
        end

        function solvesystemwithcolumn(testCase,systems,blocked,recursive2)
            options.blocked = blocked;
            options.recursive = recursive2;
            options.algorithm = 'column';
            [~,output] = macaulaylab(systems,20,options);
            result = output.accuracy;
            testCase.verifyLessThan(result,10e-12);
        end

%         function solvesystemwithnullcheb(testCase,systemscheb, ...
%                 blocked,recursive1)
%             options.blocked = blocked;
%             options.recursive = recursive1;
%             options.algorithm = 'null';
%             options.basis = 'chebyshev';
%             [~,output] = macaulaylaby(systemscheb,20,options);
%             result = output.accuracy;
%             testCase.verifyLessThan(result,10e-12);
%         end

function solvemep(testCase,l,n,dmax,blocked,recursive1)
            options.blocked = blocked;
            options.recursive = recursive1;
            [~,output] = macaulaylab(randommep(dmax,n,l+n-1,l),30,options);
            result = min(output.residuals);
            testCase.verifyLessThan(result,10e-7);
        end
    end

end