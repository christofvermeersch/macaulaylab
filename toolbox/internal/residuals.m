function [maxres, res] = residuals(solutions, problem, varargin)  
    %RESIDUALS - Mixed residuals for a problem.
    %   maxres = RESIDUALS(problem, solution) computes the maximum mixed 
    %   residual of the computed solution for a given problem.
    %
    %   [..., res] = RESIDUALS(...) also returns the individual residual
    %   value of every solution.
    %
    %   [...] = RESIDUALS(..., basis) uses a user-specified monomial basis.
    %   The default basis is the standard monomial basis.
    %
    %   Input arguments:
    %       - solutions (double): matrix of solutions.
    %       - problem (problemstruct or subclass): polynomial problem.
    %       - type (string = "mixed" - optional): residual type 
    %           ("absolute" | "relative" | "mixed").
    %       - basis (optional): monomial basis function handle.
    %
    %   Output arguments:
    %       - maxres (double): maximum residual value.
    %       - res (double): vector of residuals for each solution.
    %
    %   See also VALUE.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: added different types of residuals.
    %   - 2026 by CV: updated documentation and comments.

    % Process the optional arguments:
    type = "mixed";
    basis = @monomial;
    for i = 1:length(varargin)
        switch class(varargin{i})
            case "string"
                type = varargin{i};     
            case "function_handle"
                basis = varargin{i};
            otherwise 
                error("Argument type is not recognized.")
        end
    end

    % Compute the absolute residuals for all the solutions:
    nSolutions = size(solutions, 1);
    resabs = zeros(nSolutions, 1);
    for i = 1:nSolutions
        valuetensor = value(solutions(i, :), problem, basis);
        resabslist = zeros(problem.s, 1);
        for j = 1:problem.s
            singularvalues = svd(squeeze(valuetensor(j, :, :)));
            resabslist(j) = singularvalues(end); 
        end
        resabs(i) = sum(resabslist);
    end
    
    % Keep absolute residuals or change to relative/mixed residuals:
    if strcmp(type, "absolute")
        res = resabs;
    else
        num = resabs;
        den = 0;

        % Precompute coefficient norms (independent of solutions):
        coefnorms = cell(problem.s, 1);
        for j = 1:problem.s
            coeftensor = problem.coef{j};
            norms = zeros(size(coeftensor, 1), 1);
            for k = 1:size(coeftensor, 1)
                norms(k) = norm(squeeze(coeftensor(k, :, :)));
            end
            coefnorms{j} = norms;
        end

        for i = 1:nSolutions
            sol = solutions(i, :);
            for j = 1:problem.s
                suppmatrix = problem.supp{j};
                if isequal(basis, @monomial)
                    monomialvalues = prod(abs(sol).^suppmatrix, 2);
                    den = den + coefnorms{j}.'*monomialvalues;
                else
                    for k = 1:size(suppmatrix, 1)
                        den = den + abs(value(sol, coefnorms{j}(k), ...
                            suppmatrix(k, :), basis));
                    end
                end
            end
        end

        if strcmp(type, "relative")
            den = den/problem.s;
        elseif strcmp(type, "mixed")
            den = (den + 1)/problem.s;
        else
            error("Residual type is not recognized.")
        end
        res = num./den;
    end
    
    % Determine the maximum residual:
    maxres = max(res);
end
