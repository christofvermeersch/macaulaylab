function C = shift(A, B, basis)
    %SHIFT - Result of shifting two monomials.
    %   C = SHIFT(A, B) computes the shift of the monomials/polynomial in B 
    %   by the monomials in A in the standard monomial basis. 
    %
    %   C = SHIFT(..., basis) considers the shift in a user-specified
    %   monomial basis. This function serves as an interface for the
    %   different implemented basis operations.
    %
    %   Input arguments:
    %       - A (int): exponent matrix of the shift monomials.
    %       - B (int): exponent matrix of the polynomial to shift.
    %       - basis (function_handle = @monomial - optional): monomial 
    %           basis.
    %
    %   Output arguments:
    %       - C (double): coefficient-exponent matrix of the shifted 
    %           polynomial.
    %
    %   See also MONOMIAL, CHEBYSHEV, LEGENDRE.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Set the default monomial basis:
    if nargin < 3
        basis = @monomial;
    end

    % Compute the shifted monomial/polynomial:
    C = basis(A, B); % the function is an interface to the basis operation.
end
