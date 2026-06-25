function [P, aff, hit, rem] = rowrecomb(d, m, idx, poly, varargin)
    %ROWRECOMB - Row combination matrix.
    %   P = ROWRECOMB(d, m, idx, shift) creates a sparse recombination 
    %   matrix for a problem in which rows given by the indices are shifted  
    %   by a certain shift polynomial.
    %
    %   P = ROWRECOMB(..., blocksize) considers problems with eigenvectors.
    %
    %   P = ROWRECOMB(..., basis) performs the shift in a user-specified
    %   monomial basis.
    %
    %   P = ROWRECOMB(..., order) performs the shift in a user-specified
    %   monomial order.
    %
    %   Input arguments:
    %       - d (int): degree of the Macaulay matrix.
    %       - m (int): number of variables.
    %       - idx (int): row indices to shift.
    %       - poly (double): shift polynomial in coefficient-exponent 
    %           format.
    %       - blocksize (optional): block size (eigenvector length).
    %       - basis (function_handle = @monomial - optional): monomial 
    %           basis.
    %       - order (function_handle = @grevlex - optional): monomial 
    %           order.
    %
    %   Output arguments:
    %       - P (double): sparse recombination matrix.
    %       - aff (int): number of affected rows.
    %       - hit (int): number of hit rows.
    %       - rem (int): number of remaining rows.
    %
    %   See also ROWSHIFT.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: pre-allocated construction of B2 and rep.
    %   - 2026 by CV: updated documentation and comments.

    % Process the optional arguments:
    blocksize = 1;
    basis = @monomial;
    order = @grevlex;
    for i = 1:nargin-4
        switch class(varargin{i})
            case "function_handle"
                if nargin(varargin{i}) < 2
                    order = varargin{i};
                else
                    basis = varargin{i};
                end
            case "double"
                blocksize = varargin{i};
            otherwise 
                error("Argument type is not recognized.")
        end
    end

    % Create empty recombination matrix:
    p = nbmonomials(d, m)*blocksize;
    P = zeros(p, p);

    % Create B matrix:
    P(1:length(idx), idx) = eye(length(idx));

    Sg = rowshift(d, m, idx, poly, blocksize, basis, order);
    Sgmasked = Sg; 
    Sgmasked(:, idx) = 0;
    B2 = find(any(Sgmasked, 2))';

    % Create C matrix:
    P(length(idx)+1:length(idx)+length(B2), :) = Sg(B2, :);
    rep = [idx zeros(1, length(B2))];
    repcount = length(idx);
    for k = 1:length(B2)
        rowsinsg = find(Sg(B2(k), :));
        for l = 1:length(rowsinsg)
            if ~ismember(rowsinsg(l), rep(1:repcount))
                repcount = repcount + 1;
                rep(repcount) = rowsinsg(l);
                break
            end
        end
    end
    notrep = setdiff(1:p, rep(1:repcount));

    % Create D matrix:
    P(length(idx)+length(B2)+1:end, notrep) = eye(length(notrep));  
    P = sparse(P);
    
    % Output additional information:
    aff = length(idx);
    hit = length(B2);
    rem = p-aff-hit;
end
