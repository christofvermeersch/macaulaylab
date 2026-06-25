function pos = grvlex(A)
    %GRVLEX - Graded inverse lexicographic monomial order.
    %   pos = GRVLEX(A) returns the position of a monomial in the monomial
    %   basis using the graded inverse lexicographic order.
    %
    %   Input arguments:
    %       - A (int): matrix of monomial exponent vectors with every row
    %       corresponding to a different monomial.
    %
    %   Output arguments:
    %       - pos (int): position vector in grvlex order.
    %
    %   See also POSITION.
            
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: vectorized the computation in one direction.
    %   - 2026 by CV: updated documentation and comments.

    [nr, nc] = size(A);
    d = sum(A, 2);
    pos = nbmonomials(d, nc);
    P = [zeros(nr, 1), cumsum(A, 2)];
    for l = nc:-1:1
        j = nc - l;
        pos = pos - nbmonomials(P(:, j+1) - 1, j);
    end
end