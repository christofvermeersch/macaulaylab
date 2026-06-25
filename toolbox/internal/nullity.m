function n = nullity(M, tol)
    %NULLITY - Nullity of a matrix.
    %   n = NULLITY(M) returns the nullity of the matrix M. 
    %
    %   n = NULLITY(..., tol) uses a tolerance for the rank determination.
    %
    %   Input arguments:
    %       - M (double): matrix.
    %       - tol (double = 1e-10 - optional): numerical tolerance for the 
    %           rank determination.
    %
    %   Output arguments:
    %       - n (int): nullity of M.
    %
    %   See also RANK.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Determine the rank
    if nargin < 2
        r = rank(M);
    else
        r = rank(M, tol);
    end
    
    % Determine the nullity:
    n = size(M, 2) - r;
end