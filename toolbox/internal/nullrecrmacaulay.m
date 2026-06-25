function Z = nullrecrmacaulay(Y, N, tol)
    %NULLRECRMACAULAY - Recursive null space update of a Macaulay matrix.
    %   Z = NULLRECRMACAULAY(Y ,N) builds a numerical basis of the 
    %   null space  of the Macaulay  matrix, using a  numerical basis of 
    %   the null space of a lower degree Macaulay matrix and the difference 
    %   between both Macaulay matrices.
    %
    %   Z = NULLRECRMACAULAY(..., tol) uses a user-specified tolerance for
    %   the rank decision.
    %
    %   Input arguments:
    %       - Y (double): basis matrix of the null space at the lower 
    %           degree.
    %       - N (double): new rows added to the Macaulay matrix.
    %       - tol (double = 1e-10 - optional): numerical tolerance for 
    %           near-zero matrices.
    %
    %   Output arguments:
    %       - Z (double): updated basis matrix of the null space.
    %
    %   See also NULL, MACAULAYUPDATE, NULLSPARSEMACAULAY,
    %   NULLITERMACAULAY.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.
    
    % Process the optional parameter:
    if nargin < 3
        tol = 1e-10;
    end

    % Compute the updated basis matrix:
    [q, n] = size(Y); 
    X = N(:, 1:q)*Y;
    [~, s, V] = svd(full([X N(:, q+1:end)]));
    r = nnz(s > tol);
    V = V(:, r+1:end);
    Z = [Y*V(1:n, :); V(n+1:end, :)];
end

