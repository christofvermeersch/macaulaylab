function pos = grnlex(A)
    %GRNLEX - Graded negative lexicographic monomial order.
    %   pos = GRNLEX(A) returns the position of a monomial in the monomial
    %   basis using the graded negative lexicographic order. 
    %
    %   Input arguments:
    %       - A (int): matrix of monomial exponent vectors with every row
    %       corresponding to a different monomial.
    %
    %   Output arguments:
    %       - pos (int): position vector in grnlex order.
    %
    %   See also POSITION.
            
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: vectorized the computation in one direction.
    %   - 2026 by CV: updated documentation and comments.

    [nr, nc] = size(A);
    S = fliplr(cumsum(fliplr(A), 2));
    pos = zeros(nr, 1);
    for l = 1:nc-1
        pos = pos + nbmonomials(S(:, l) - 1, nc - l + 1);
    end
    pos = pos + A(:, nc) + 1;
end