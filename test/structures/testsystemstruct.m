classdef testsystemstruct < matlab.unittest.TestCase
    methods (Test)
        % Test the class constructor:
        function constructortest(testCase)
            p = [1 5 0; 2 2 2];
            q = [3 4 2; 4 1 1];
            system = systemstruct({p,q});
            testCase.verifyClass(system,"systemstruct");
            testCase.verifyInstanceOf(system,"problemstruct");
            testCase.verifyEqual(system.eqs,{p,q});
            testCase.verifyEqual(system.coef,{[1; 2];[3; 4]});
            testCase.verifyEqual(system.supp,{[5 0; 2 2];[4 2; 1 1]});
            testCase.verifyEqual(system.basis,NaN);
        end

        function symbolictest(testCase)
            syms x y
            p = x^5 + 2*x^2*y^2;
            q = 3*x^4*y^2 + 4*x*y;
            system = systemstruct({p, q});
            testCase.verifyClass(system, "systemstruct");
            testCase.verifyInstanceOf(system, "problemstruct");
            testCase.verifyEqual(system.coef, {[1; 2]; [3; 4]});
            testCase.verifyEqual(system.supp, {[5 0; 2 2]; [4 2; 1 1]});
        end

        function maptest(testCase)
            p = [1 5 0; 2 2 2];
            q = [3 4 2; 4 1 1];
            system = systemstruct({p, q});
            testCase.verifyEqual(system.map, NaN);

            syms x y
            p2 = x^5 + 2*x^2*y^2;
            q2 = 3*x^4*y^2 + 4*x*y;
            systemsym = systemstruct({p2, q2});
            testCase.verifyEqual(systemsym.map, [x y]);
        end
    end
end
