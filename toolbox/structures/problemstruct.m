classdef problemstruct
    %PROBLEMSTRUCT - Class for generic problems.
    %   obj = PROBLEMSTRUCT(coef, supp) creates a structure that contains
    %   all the information of the problem defined by the coefficients and
    %   support. 
    %
    %   obj = PROBLEMSTRUCT(..., basis) also specifies the monomial basis
    %   of the problem.
    %
    %   Notice that every matrix polynomial has one entry in the s-by-1 
    %   cells coef and supp. The coefficients or k-by-l coefficient 
    %   matrices are given in a r-by-k-by-l tensor and every mode 
    %   corresponds to one row in the r-by-n support matrix.
    %
    %   List of properties:
    %       - coef (double): tensor with the coefficient matrices.
    %       - supp (int): support of the coefficient matrices.
    %       - basis (function_handle - optional): monomial basis.
    %       - m (int): number of variables/parameters.
    %       - s (int): number of (matrix) equations.
    %       - k (int): number of rows of the coefficient matrices.
    %       - l (int): number of columns of the coefficient matrices.
    %       - di (cell): total degree of every equation.
    %       - dmax (int): maximum total degree in the support.
    %       - nnze (cell): number of non-zero coefficients.
    %       - ma (int): number of affine solutions.
    %       - mb (int): total number of solutions.
    %       - multiplicities (logical): (logical): flag for multiplicities.
    %       - posdim (logical): flag for positive-dimensional solution 
    %           sets.
    %       - solved (logical): flag set to indicatate whether the problem 
    %           is solved in the database.
    %       - appid (string): related applications.
    %       - bibid (string): references in BibTeX.
    %       - comments (string): additional comments about the problem.
    %
    %   List of methods:
    %       - properties: returns the most important properties.
    %       - disp: displays a short description.
    %       - info: gives a more elaborate description.
    %       - structure: visualizes the problem structure (support).
    %       - coefmat: creates the coefficient matrices.
    %       - coefeqs: return the polynomial coefficients.
    %
    %   See also SYSTEMSTRUCT, MEPSTRUCT.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.
    %   - 2026 by CV: changed number of parameters to variable "m".
    %   - 2026 by CV: added coefeqs method to function.

    properties
        coef
        supp
        m
        s
        k
        l
        di
        dmax
        nnze
        basis
        ma
        mb
        multi
        posdim
        solved
        appid
        bibid
        comments
    end

    methods
        function obj = problemstruct(coef, supp, basis)
            %PROBLEMSTRUCT - Constructor for the problemstruct class.
            %   obj = PROBLEMSTRUCT(coef, supp) creates a structure that
            %   contains all the information of the problem defined by the
            %   coefficients and support.
            %
            %   obj = PROBLEMSTRUCT(..., basis) also specifies the monomial
            %   basis of the problem.
            %
            %   Input arguments:
            %       - coef (cell): cell array with coefficient tensors.
            %       - supp (cell): cell array with support matrices.
            %       - basis (function_handle - optional): monomial basis.
            %
            %   Output arguments:
            %       - obj (problemstruct): problem structure.

            % Process the optional argument:
            if nargin < 3
                obj.basis = NaN;
            else
                obj.basis = basis;
            end
            
            % Create a problemstruct:
            obj.coef = coef;
            obj.supp = supp;
            obj.m = size(supp{1}, 2);
            obj.s = numel(coef);
            obj.di = cell(obj.s, 1);
            obj.nnze = cell(obj.s, 1);
            for k = 1:obj.s
                obj.di{k} = max(sum(supp{k}, 2));
                obj.nnze{k} = nnz(coef{k});
            end   
            obj.dmax = max(horzcat(obj.di{:}));
            [~, obj.k, obj.l] = size(coef{1});
        end 

        function [s, dmax, m, di, nnze, k, l] = properties(obj)
            %PROPERTIES - Most important properties of the problem.
            %   [...] = PROPERTIES(problem) returns the
            %   most important properties of the problem.
            %
            %   Input arguments:
            %       - obj (problemstruct): polynomial problem structure.
            %
            %   Output arguments:
            %       - s (int): number of (matrix) equations.
            %       - dmax (int): maximum total degree.
            %       - m (int): number of variables.
            %       - di (cell): total degrees per equation.
            %       - nnze (cell): non-zero counts per equation.
            %       - k (int): number of rows of the coefficient matrices.
            %       - l (int): number of columns of the coefficient matrices.

            s = obj.s;
            m = obj.m;
            di = obj.di;
            dmax = obj.dmax;
            nnze = obj.nnze;
            k = size(obj.coef{1}, 2);
            l = size(obj.coef{1}, 3);
        end

        function output = disp(obj)
            %DISP - Description of the problem.
            %   DISP(problem) displays a short description of the problem.
            %
            %   Input arguments:
            %       - obj (problemstruct): polynomial problem structure.
            %
            %   Output arguments:
            %       - output (string): description (when nargout > 0).
        
            switch class(obj)
                case "systemstruct"
                    if obj.s < 2
                        if obj.m < 2
                            text = sprintf("     univariate polynomial" ...
                                + " with degree %d \n", obj.dmax);
                        else
                            text = sprintf("     multivariate " ...
                                + "polynomial with total degree %d \n", ...
                                obj.dmax);
                        end
                    else 
                       if obj.m < 2
                            text = sprintf("     system of %d " ... 
                                + "univariate polynomials with maximum" ...
                                + " degree %d \n", obj.s, obj.dmax);
                        else
                            text = sprintf("     system of %d " ... 
                                + "multivariate polynomials with " ...
                                + "maximum total degree %d \n", obj.s, ...
                                obj.dmax);
                        end
                    end
                    
                case "mepstruct"
                    if obj.m < 2
                        text = sprintf("     %d-parameter eigenvalue " ...
                            + "problem with degree %d \n", obj.m, ...
                            obj.dmax);
                    else
                        text = sprintf("     %d-parameter eigenvalue " ...
                            + "problem with total degree %d \n", ...
                            obj.m, obj.dmax);
                    end
                otherwise
                    if obj.k > 1 || obj.l > 1
                        text = sprintf("     problem with %d matrix " ...
                            + "equations in %d variables \n", obj.s, obj.m);
                    else 
                        text = sprintf("     problem with %d equations" ...
                            + " in %d variables \n", obj.s, obj.m);
                    end
            end
            
            if nargout > 0
                output = text;
            else
                fprintf(text);
            end
        end

        function output = info(obj)
            %INFO - Information about the problem.
            %   INFO(problem) displays a more elaborate description of the
            %   problem that contains all the relevant information.
            %
            %   Input arguments:
            %       - obj (problemstruct): problem structure.
            %
            %   Output arguments:
            %       - output (string): description (when nargout > 0).

            text = disp(obj);
            text = text + structure(obj);
            text = text + sprintf("    sparsity factor is %.2f \n", sparsity(obj));

            if nargout > 0
                output = text;
            else
                fprintf(text);
            end
        end

        function output = structure(obj)
            %STRUCTURE - Visualization of the problem structure.
            %   STRUCTURE(problem) visualizes the problem structure.
            %
            %   Input arguments:
            %       - obj (problemstruct): polynomial problem structure.
            %
            %   Output arguments:
            %       - output (string): visualization (when nargout > 0).

            if obj.k > 1 || obj.l > 1
                text = newline;
                for i = 1:obj.s
                    support = obj.supp{i};
                    text = text + sprintf("        (");
                    for j = 1:size(support, 1)
                        text = text + sprintf("A");
                        for vi = 1:obj.m
                            text = text + sprintf("%d", support(j, vi));
                        end
                        text = text + sprintf(" ");
                        if sum(support(j, :)) > 0
                            for vi = 1:obj.m
                                if support(j, vi) > 0
                                    text = text + sprintf("l%d", vi);
                                end
                                if support(j, vi) > 1
                                    text = text + sprintf("^%d ", ...
                                        support(j, vi));
                                end
                            end
                        end
                        if j < size(support, 1)
                            text = text + sprintf("+ ");
                        end
                    end
                    text = text + sprintf(") z = 0 \n");
                end
                text = text + newline;
                text = text + sprintf("    coefficient matrices are " ...
                    + "%d x %d matrices \n", obj.k, obj.l);
            else 
                text = newline;
                for i = 1:obj.s
                    text = text + sprintf("        p%d(x1, ..., x%d) = 0" ...
                        + "     <- d%d = %d \n", i, obj.m, i, obj.di{i});
                end
                text = text + newline;
            end

            if nargout > 0
                output = text;
            else
                fprintf(text);
            end
        end

        function mat = coefmat(obj)
            %COEFMAT - Coefficient matrices of the matrix equation.
            %   mat = COEFMAT(problem) creates the coefficient matrices
            %   from the coefficient tensor (only works for single matrix
            %   equations).
            %
            %   Input arguments:
            %       - obj (problemstruct): problem structure for a problem
            %           with a single matrix equation.
            %
            %   Output arguments:
            %       - mat (cell): coefficient matrices.

            if obj.s > 1
                error("Method only works for single matrix equations.")
            end
            nMat = size(obj.supp{1}, 1);
            mat = cell(nMat, 1);
            tens = obj.coef{1};
            for i = 1:nMat
                mat{i} = squeeze(tens(i, :, :));
            end
        end

        function eqs = coefeqs(obj)
            %COEFEQS - Coefficients of the multivariate polynomial system.
            %   eqs = COEFEQS(problem) creates the coefficients from the
            %   coefficient tensor (only works for multivariate polynomial
            %   systems).
            %
            %   Input arguments:
            %       - obj (problemstruct): problem structure for a problem 
            %           with scalar polynomials.
            %
            %   Output arguments:
            %       - eqs (cell): coefficient vectors.

            if obj.k > 1 || obj.l > 1
                error("Method only works for scalar polynomials")
            end
            eqs = cell(obj.s, 1);
            for i = 1:obj.s
                eqs{i} = squeeze(obj.coef{i});
            end
        end
    end
end
