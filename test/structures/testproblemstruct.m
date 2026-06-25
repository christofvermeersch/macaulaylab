classdef testproblemstruct < matlab.unittest.TestCase
    methods (Test)
        % Test the class constructor:
        function constructortest(testCase)
            coef = {randn(4,3,2);randn(4,3,2)};
            supp = {[0 0; 1 0; 0 1; 2 0];[0 0; 2 2; 0 1; 0 5]};
            problem = problemstruct(coef,supp);
            testCase.verifyClass(problem,"problemstruct");
            testCase.verifyEqual(problem.coef,coef);
            testCase.verifyEqual(problem.supp,supp);
            testCase.verifyEqual(problem.m,2);
            testCase.verifyEqual(problem.s,2);
            testCase.verifyEqual(problem.k,3);
            testCase.verifyEqual(problem.l,2);
            testCase.verifyEqual(problem.di,{2;5});
            testCase.verifyEqual(problem.dmax,5);
            testCase.verifyEqual(problem.nnze,{24;24});
            testCase.verifyEqual(problem.basis,NaN);
        end

        % Test the class functions:
        function propertiestest(testCase)
            coef = {[1; 2];[3; 4]};
            supp = {[5 0; 2 2];[4 2; 1 1]};
            problem = problemstruct(coef,supp);
            [s,dmax,m,di,nnze,k,l] = properties(problem);
            expected = {2,6,2,{5;6},{2;2},1,1};
            testCase.verifyEqual({s,dmax,m,di,nnze,k,l},expected);
        end

        function constructorbasistest(testCase)
            coef = {[1; 2]};
            supp = {[3 0; 1 1]};
            problem = problemstruct(coef, supp, @chebyshev);
            testCase.verifyEqual(problem.basis, @chebyshev);
        end

        function dispsystemtest(testCase)
            p = [1 3; 2 0];
            out = disp(systemstruct({p}));
            testCase.verifyTrue(contains(out, "univariate polynomial"));

            p = [1 1 1; 2 2 0];
            out = disp(systemstruct({p}));
            testCase.verifyTrue(contains(out, "multivariate polynomial"));

            p = [1 3; 2 0]; q = [1 1; 3 0];
            out = disp(systemstruct({p, q}));
            testCase.verifyTrue(contains(out, "univariate polynomials"));

            p = [1 1 1; 2 2 0]; q = [3 1 0; 1 0 1];
            out = disp(systemstruct({p, q}));
            testCase.verifyTrue(contains(out, "multivariate polynomials"));
        end

        function dispmeptest(testCase)
            A0 = [1 0; 0 1]; A1 = [0 1; 1 0];
            out = disp(mepstruct({A0, A1}, [0; 1]));
            testCase.verifyTrue(contains(out, "parameter eigenvalue problem"));

            out = disp(mepstruct({A0, A1, A0}, [0 0; 0 1; 1 0]));
            testCase.verifyTrue(contains(out, "parameter eigenvalue problem"));
        end

        function dispproblemstructtest(testCase)
            coef = {[1; 2]; [3; 4]};
            supp = {[5 0; 2 2]; [4 2; 1 1]};
            out = disp(problemstruct(coef, supp));
            testCase.verifyTrue(contains(out, "equations"));

            coef = {randn(2, 3, 2)};
            supp = {[0 0; 1 0]};
            out = disp(problemstruct(coef, supp));
            testCase.verifyTrue(contains(out, "matrix equations"));
        end

        function infotest(testCase)
            p = [1 1 1; 2 2 0]; q = [3 1 0; 1 0 1];
            system = systemstruct({p, q});
            out = info(system);
            testCase.verifyClass(out, "string");
        end

        function structuretest(testCase)
            p = [1 1 1; 2 2 0]; q = [3 1 0; 1 0 1];
            system = systemstruct({p, q});
            out = structure(system);
            testCase.verifyClass(out, "string");
            testCase.verifyTrue(contains(out, "p1"));

            coef = {randn(2, 3, 2)};
            supp = {[0 0; 1 0]};
            problem = problemstruct(coef, supp);
            out = structure(problem);
            testCase.verifyClass(out, "string");
            testCase.verifyTrue(contains(out, "A"));
        end

        function coefmattest(testCase)
            A0 = [1 2; 3 4]; A1 = [0 1; 1 0];
            mep = mepstruct({A0, A1}, [0; 1]);
            mat = coefmat(mep);
            testCase.verifyClass(mat, "cell");
            testCase.verifyEqual(mat{1}, A0);
            testCase.verifyEqual(mat{2}, A1);

            coef = {randn(2, 3, 2); randn(2, 3, 2)};
            supp = {[0 0; 1 0]; [0 0; 0 1]};
            problem = problemstruct(coef, supp);
            testCase.verifyError(@() coefmat(problem), '');
        end

        function coefeqstest(testCase)
            p = [1 5 0; 2 2 2]; q = [3 4 2; 4 1 1];
            system = systemstruct({p, q});
            eqs = coefeqs(system);
            testCase.verifyClass(eqs, "cell");
            testCase.verifyEqual(eqs{1}, [1; 2]);
            testCase.verifyEqual(eqs{2}, [3; 4]);

            coef = {randn(2, 3, 2)};
            supp = {[0 0; 1 0]};
            problem = problemstruct(coef, supp);
            testCase.verifyError(@() coefeqs(problem), '');
        end
    end
end