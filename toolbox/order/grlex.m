function pos = grlex(A)
    %GRLEX - Graded lexicographic monomial order.
    %   pos = GRLEX(A) returns the position of a monomial in the monomial
    %   basis using the graded lexicographic order. 
    %
    %   Input arguments:
    %       - A (int): matrix of monomial exponent vectors with every row
    %       corresponding to a different monomial.
    %
    %   Output arguments:
    %       - pos (int): position vector in grlex order.
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
    S = [fliplr(cumsum(fliplr(A), 2)), zeros(nr, 1)];
    for l = 1:nc
        pos = pos - nbmonomials(S(:, l+1) - 1, nc - l);
    end
end