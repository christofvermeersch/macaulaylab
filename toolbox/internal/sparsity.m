function s = sparsity(problem, zerotol)
    %SPARSITY - Sparsity factor of a problem.
    %   s = SPARSITY(problem) returns the sparsity factor of a problem,
    %   which is the number of non-zero elements divided by the possible
    %   number of elements if all coefficient (matrices) would be non-zero.
    %   
    %   s = SPARSITY(..., zerotol) ignores all coefficients smaller than 
    %   the given tolerance.
    %
    %   Input arguments:
    %       - problem (problemstruct or subclass): polynomial problem.
    %       - zerotol (double = 1e-16 - optional): tolerance for
    %           considering coefficients equal to zero.
    %
    %   Output arguments:
    %       - s (double): sparsity factor (non-zero count divided by 
    %           maximum possible count).

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Process the optional argument:
    if nargin < 2
        zerotol = 1e-16;
    end

    % Compute the sparsity factor:
    num = 0;
    den = 0;
    for i = 1:problem.s
        nnze = nnz(abs(problem.coef{i}) >= zerotol);
        num = num + nnze;
        den = den + problem.k*problem.l* ...
            nbmonomials(problem.di{i}, problem.m);
    end
    s = num/den;
end
