function pos = grevlex(A)
    %GREVLEX   Graded reverse lexicographic monomial order.
    %   pos = GREVLEX(A) returns the position of a monomial in the monomial 
    %   basis using the graded reverse lexicographic order. The result
    %   is a k-by-1 vector that contains the position for every n-variate
    %   monomial in the k-by-n matrix A. 
    %
    %   See also POSITION.
            
    % Copyright (c) 2025 - Christof Vermeersch
    % 
    % Updates:
    %   - 2025: vectorized the computation in one direction.
    
    [nr,nc] = size(A);
    B = cumsum(A,2)-1;
    pos = ones(nr,1);
    for l = nc:-1:1
        pos = pos + nbmonomials(B(:,nc-l+1),nc-l+1);
    end
end