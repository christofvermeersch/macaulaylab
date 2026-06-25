function [C, alpha, beta, phi0, phi1] = monomial(A, B)
    %MONOMIAL - Standard multivariate monomial basis.
    %   C = MONOMIAL(A, B) computes the product of the two monomials. If A 
    %   and/or B are polynomials, then all monomials of B are multiplied by
    %   every monomial of A and stacked on top of each other in the 
    %   polynomial C.
    %
    %   [..., alpha, beta, phi0, phi1] = MONOMIAL(...) also returns the 
    %   alpha(x) and beta(x) of the recursion relation for the standard 
    %   monomial basis, i.e., x^{k + 1} = alpha(x) * x^k + beta(x) 
    %   x^{k - 1}. The first two basis functions are phi0(x) and phi1(x).
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

    % Perform the multiplication:
    C = [ones(size(A, 1)*size(B, 1), 1) ...
        kron(A, ones(size(B, 1), 1)) + kron(ones(size(A, 1), 1), B)];

    % Set the properties of the monomial basis:
    n = size(A, 2);
    alpha = ones(1, n+1);
    beta = [0 zeros(1, n)];
    phi0 = [1 zeros(1, n)];
    phi1 = ones(1, n+1);
end
