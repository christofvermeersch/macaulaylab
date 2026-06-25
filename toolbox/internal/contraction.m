function q = contraction(p, d, m, varargin)
    %CONTRACTION - Matrix representation of a polynomial.
    %   q = CONTRACTION(p, d, m) contracts a polynomial given by its vector 
    %   representation (in the graded reverse lexicographic ordering) into 
    %   its corresponding coefficient matrix.
    %
    %   q = CONTRACTION(..., tol) only considers coefficients larger
    %   than the user-specified tolerance.
    %
    %   q = CONTRACTION(..., order) uses a user-specified monomial order.
    %
    %   Input arguments:
    %       - p (double): polynomial in vector representation.
    %       - d (int): maximum total degree.
    %       - m (int): number of variables.
    %       - tol (double = 1e-10 - optional): numerical tolerance for zero 
    %           coefficients.
    %       - order (function_handle = @grevlex - optional): monomial 
    %           order.
    %
    %   Output arguments:
    %       - q (double): polynomial in coefficient-exponent matrix 
    %           representation.
    %
    %   See also EXPANSION.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Process the optional parameters:
    tol = 1e-10;
    order = @grevlex;
    for i = 1:nargin-3
        switch class(varargin{i})
            case "function_handle"
                order = varargin{i};
            case "double"
                tol = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end
       
    % Contract the polynomial:
    p = fliplr(p);
    mask = find(abs(p) > tol);
    K = monomials(d, m, order);
    q = [p(mask).' K(mask, :)];    
end
