function s = nbmonomials(d, m, blocksize)
    %NBMONOMIALS - Number of monomials.
    %   s = NBMONOMIALS(d, m) returns the number of monomials in the
    %   monomial basis.
    %
    %   s = NBMONOMIALS(..., blocksize) takes into account the blocksize.
    %
    %   Input arguments:
    %       - d (int): maximum total degree.
    %       - m (int): number of variables.
    %       - blocksize (int = 1 - optional): block size.
    %
    %   Output arguments:
    %       - s (int): number of monomials.
    %
    %   See also MONOMIALS.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.
    
    % Compute the number of monomials:
    if d < 0
        s = 0; % avoid errors with negative degrees
    else
        s = round(exp(gammaln(d+m+1)-gammaln(m+1)-gammaln(d+1)));
    end

    % Take into account the blocksize:
    if nargin > 2
        s = s*blocksize;
    end
end
