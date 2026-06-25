function pos = grevlex(A)
    %GREVLEX - Graded reverse lexicographic monomial order.
    %   pos = GREVLEX(A) returns the position of a monomial in the monomial
    %   basis using the graded reverse lexicographic order.
    %
    %   Input arguments:
    %       - A (int): matrix of monomial exponent vectors with every row
    %       corresponding to a different monomial.
    %
    %   Output arguments:
    %       - pos (int): position vector in grevlex order.
    %
    %   See also POSITION.
            
    % Copyright (c) 2026 - Christof Vermeersch
    % 
    % Updates:
    %   - 2025 by CV: vectorized the computation in one direction.
    %   - 2026 by CV: updated documentation and comments.
    
    [nr, nc] = size(A);
    B = cumsum(A, 2)-1;
    pos = ones(nr, 1);
    for l = nc:-1:1
        pos = pos + nbmonomials(B(:, nc-l+1), nc-l+1);
    end
end