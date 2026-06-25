function [C, alpha, beta, phi0, phi1] = legendre(A, B)
    %LEGENDRE - Multivariate product Legendre basis.
    %   C = LEGENDRE(A, B) computes the product of the two monomials. If A 
    %   and/or B are polynomials, then all monomials of B are multiplied by
    %   every monomial of A and stacked on top of each other in the 
    %   polynomial C.
    %
    %   [..., alpha, beta, phi0, phi1] = LEGENDRE(...) also returns the  
    %   alpha(x) and beta(x) of the recursion relation for the Legendre 
    %   basis, i.e., P_{k + 1}(x) = alpha(x) * P_k(x) + beta(x) 
    %   P_{k - 1}(x). The first two basis functions are given in phi0(x)
    %   and phi1(x).
    %
    %   Note that this function is not yet implemented.
    %
    %   Input arguments:
    %       - A (int): exponent matrix of the first polynomial.
    %       - B (int): exponent matrix of the second polynomial.
    %
    %   Output arguments:
    %       - C (double): coefficient-exponent matrix of the product.
    %       - alpha (double): multiplication coefficient vector of the 
    %           recursion.
    %       - beta (double): addition coefficient vector of the recursion.
    %       - phi0 (double): first basis function.
    %       - phi1 (double): second basis function.
    %
    %   See also SHIFT.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    error("Basis operation is not yet implemented.")    
end    
