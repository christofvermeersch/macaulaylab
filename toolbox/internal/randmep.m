function mep = randmep(dmax, m, k, l, varargin)
    %RANDMEP - Dense multiparameter eigenvalue problem.
    %   mep = RANDMEP(dmax, m, k, l) generates a dense multiparameter
    %   eigenvalue problem.
    %
    %   mep = RANDMEP(..., complex) adds the option for the coefficient
    %   matrices to be complex (= boolean).
    %
    %   mep = RANDMEP(..., supp) results in a multiparameter eigenvalue
    %   problem with only non-zero coefficient matrices for the elements in
    %   the given support.
    %
    %   Input arguments:
    %       - dmax (int): maximum total degree.
    %       - m (int): number of eigenvalue parameters.
    %       - k (int): number of rows of the coefficient matrices.
    %       - l (int): number of columns of the coefficient matrices.
    %       - complex (logical = false - optional): flag for complex 
    %           coefficient matrices.
    %       - supp (cell - optional): support matrix.
    %
    %   Output arguments:
    %       - mep (mepstruct): multiparameter eigenvalue problem.
    %
    %   See also MEPSTRUCT.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Process the optional parameters:
    basis = "unknown";
    isComplex = false;
    supp = monomials(dmax, m);
    for i = 1:nargin-4
        switch class(varargin{i})
            case "function_handle"
                basis = varargin{i};
            case "logical"
                isComplex = varargin{i};
            case "double"
                supp = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end

    % Create the mepstruct:
    dmax = max(sum(supp, 2));
    s = nbmonomials(dmax, m);
    mat = cell(s, 1);
    mat(:) = {zeros(k, l)};
    if ~isComplex % for real coefficient matrices.
        for i = grevlex(supp)'
            mat{i} = randn(k, l);
        end
    else % for complex coefficient matrices.
        for i = grevlex(supp)'
            mat{i} = randn(k, l) + 1i*randn(k, l);
        end
    end
    mep = mepstruct(mat, supp, basis);
end