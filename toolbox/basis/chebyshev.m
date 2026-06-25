function [C, alpha, beta, phi0, phi1] = chebyshev(A, B)
    %CHEBYSHEV - Multivariate product Chebyshev basis.
    %   C = CHEBYSHEV(A,B) computes the product of the two monomials. If A 
    %   and/or B are polynomials, then all monomials of B are multiplied by
    %   every monomial of A and stacked on top of each other in the 
    %   polynomial C.
    %
    %   [..., alpha, beta, phi0, phi1] = CHEBYSHEV(...) also returns the 
    %   alpha(x) and beta(x) of the recursion relation for the Chebyshev 
    %   basis, i.e., T_{k + 1}(x) = alpha(x) * T_k(x) + beta(x) 
    %   T_{k - 1}(x). The first two basis functions are given in phi0(x)
    %   and phi1(x).
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
    n = size(A, 2);
    power = 2^n;
    nMonomials = size(A, 1)*size(B, 1);
    C = zeros(power*nMonomials, n+1);
    for i = 1:size(A, 1) % loops over different monomials in A.
        D = B;
        for j = 1:n % loops over the shift in every variable.
            aij = [zeros(1, j-1) A(i, j) zeros(1, n-j)];
            D = [aij + D; abs(aij - D)];
        end
        C((i-1)*power*size(B, 1)+1:i*power*size(B, 1), 2:end) = D;
    end
    C(:, 1) = 1/power;

    % Set the properties of the monomial basis:
    alpha = [2 ones(1, n)];
    beta = [-1 zeros(1, n)];
    phi0 = [1 zeros(1, n)]; 
    phi1 = ones(1, n+1);
end
