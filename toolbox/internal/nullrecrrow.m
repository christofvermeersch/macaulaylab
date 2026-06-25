function Z = nullrecrrow(Y, N, tol)
    %NULLRECRROW - Recursive null space update of a (block) row matrix.
    %   Z = NULLRECRROW(Y, N) builds a numerical basis of the null space 
    %   of a (block) row  matrix of degree, using a numerical basis of 
    %   the null space of a lower degree (block) row matrix and the 
    %   difference between both (block) row matrices.
    %
    %   Z = NULLRECRROW(..., tol) uses a user-specified tolerance to catch
    %   near-zero matrices.
    %
    %   Input arguments:
    %       - Y (double): basis matrix of the null space at the lower 
    %           degree.
    %       - N (double): new (block) row added to the matrix.
    %       - tol (double = 1e-10 - optional): numerical tolerance for 
    %           near-zero matrices.
    %
    %   Output arguments:
    %       - Z (double): updated basis matrix of the null space.
    %
    %   See also NULL.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.
    
    % Process the optional parameter:
    if nargin < 3
        tol = 1e-10;
    end

    % Compute the updated basis matrix:
    X = N*Y;
    [~, S, V] = svd(X);
    r = nnz(diag(S) > tol);
    if r == 0
        V = eye(size(V, 1));
    end
    Z = Y*V(:, r+1:end);
end