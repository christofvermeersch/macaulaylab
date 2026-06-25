function Z = columncompr(Z, dcompr, m, varargin)
    %COLUMNCOMPR - Column compression.
    %   Z = COLUMNCOMPR(Z, dcompr, m) compresses a basis matrix of the
    %   (right) null space of a (block) Macaulay matrix with blocksize 
    %   equal to one.
    %
    %   Z = COLUMNCOMPR(..., blocksize) considers a given block size.
    %
    %   Z = COLUMNCOMPR(..., tol) uses a user-defined tolerance.
    %
    %   Input arguments:
    %       - Z (double): basis matrix of the null space to compress.
    %       - dcompr (int): compression degree.
    %       - m (int): number of variables.
    %       - blocksize (int = 1 - optional): size of the eigenvector.
    %       - tol (double = 1e-10 - optional): numerical tolerance for the 
    %           rank decision.
    %
    %   Output arguments:
    %       - Z (double): compressed basis matrix of the null space.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Process the optional parameters:
    blocksize = 1;
    tol = 1e-10;
    for i = 1:nargin-3
        if strcmp(class(varargin{i}), "double")  
            if varargin{i} < 1
                tol = varargin{i};
            else
                blocksize = varargin{i};
            end
        else
            error("Argument type is not recognized.")
        end
    end
    
    % Perform the column compression:
    nbrows = blocksize*nbmonomials(dcompr, m);
    [U, S, ~] = svd(Z(1:nbrows, :));
    ma = nnz(diag(S) > tol);
    Z = U(:, 1:ma);
end
