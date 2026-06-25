function [dgap, ma, gapsize] = gapstdmonomials(c, d, m, blocksize)
    %GAPSTDMONOMIALS - Gap zone via the standard monomials.
    %   dgap = GAPSTDMONOMIALS(c, d, m) returns the first degree of the gap
    %   zone in the list of standard monomials or NaN if there is no gap
    %   zone.
    %
    %   [dgap, ma, gapsize] = GAPSTDMONOMIALS(...) also returns the number
    %   of affine solutions and the size of the gap zone.
    %
    %   [...] = GAPSTDMONOMIALS(..., blocksize) considers the number of
    %   rows in the block rows (cf., the null space of block Macaulay
    %   matrix).
    %
    %   Input arguments:
    %       - c (int): vector of standard monomial positions.
    %       - d (int): maximum total degree to consider.
    %       - m (int): number of variables.
    %       - blocksize (int = 1 - optional): size of eigenvector.
    %
    %   Output arguments:
    %       - dgap (int): first degree of the gap zone (NaN if none).
    %       - ma (int): number of affine solutions.
    %       - gapsize (int): size of the gap zone (in number of degree 
    %           blocks).
    %
    %   See also STDMONOMIALS, COLUMNCOMPR, GAP.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Process the optional parameter:
    if nargin < 4
        blocksize = 1;
    end
    
    % Determine gap, number of affine solutions, and size of the gap:
    ma = 0;
    gapsize = 0;
    dgap = NaN;
    for k = 0:d
        nm = blocksize*nbmonomials(k, m);
        r = sum(c<=nm);
        if ma == r
            gapsize = gapsize + 1;
            dgap = k - gapsize + 1;
        else
            if ~isnan(dgap)
                break
            else
                ma = r;
            end
        end
    end     
end