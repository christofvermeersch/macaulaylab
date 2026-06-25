function V = vandermonde(d, m, x, varargin)
    %VANDERMONDE - Multivariate (block) Vandermonde matrix.
    %   V = VANDERMONDE(d, n, x) constructs the multivariate Vandermonde
    %   matrix in the standard monomial basis and graded reverse
    %   lexicographic order. The Vandermonde matrix is built from the
    %   given different m-variate points in x.
    %
    %   V = VANDERMONDE(..., z) considers the blocked version instead. The 
    %   elements of the vectors are also restructured as specified by the
    %   monomial order.
    %
    %   V = VANDERMONDE(..., order) works with a user-specified monomial
    %   order.
    %
    %   Input arguments:
    %       - d (int): maximum total degree.
    %       - m (int): number of variables.
    %       - x (double): valuation points (every row is a point).
    %       - z (double = ones - optional): eigenvector matrix.
    %       - order (function_handle = @grevlex - optional): monomial 
    %           order.
    %
    %   Output arguments:
    %       - V (double): (block) Vandermonde matrix.
    %
    %   See also MONOMIALSMATRIX.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Process the optional parameters:
    z = ones(1, size(x, 1)); % default is already transposed.
    order = @grevlex;
    isDifferentOrder = false;
    for i = 1:nargin-3
        switch class(varargin{i})
            case "function_handle"
                order = varargin{i};
                isDifferentOrder = true;
            case "double"
                z = varargin{i}.'; % transpose to put correctly in V.
            otherwise 
                error("Argument type is not recognized.")
        end
    end

    % Restructure the vectors:
    idx = order(eye(size(z, 1)))-1;
    z = z(idx, :);

    % Construct Vandermonde matrix:
    K = monomials(d, m);
    V = zeros(size(K, 1)*size(z, 1), size(x, 1));
    for k = 1:size(x, 1)
        V(:, k) = kron(prod(x(k, :).^K, 2), z(:, k));
    end
    
    % Change the order of the rows when a different monomial order is used:
    if isDifferentOrder
        idx = kron(order(K)-1, size(z, 1)*ones(size(z, 1), 1)) ...
            + kron(ones(1, size(K, 1)), 1:size(z, 1))';
        V = V(idx, :);
    end
end