function system = randsystem(s, dmax, m, varargin)
    %RANDSYSTEM - Dense system with random coefficients.
    %   system = RANDSYSTEM(s, dmax, m) generates a dense system of 
    %   multivariate polynomial equations.
    %
    %   system = RANDSYSTEM(..., complex) adds the option for the 
    %   coefficients to be complex.
    %
    %   system = RANDSYSTEM(..., supp) results in a multivariate polynomial
    %   system with only non-zero coefficients for the elements in the
    %   given support.
    %
    %   Input arguments:
    %       - s (int): number of equations.
    %       - dmax (int): maximum total degree.
    %       - m (int): number of variables.
    %       - complex (logical = false - optional): flag for complex 
    %           coefficients.
    %       - supp (cell - optional): support matrices per equation.
    %
    %   Output arguments:
    %       - system (systemstruct): polynomial system.
    %
    %   See also SYSTEMSTRUCT.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Process the optional parameters:
    basis = "unknown";
    isComplex = false;
    supp = cell(s, 1);
    supp(:) = {monomials(dmax, m)};
    for i = 1:nargin-3
        switch class(varargin{i})
            case "function_handle"
                basis = varargin{i};
            case "logical"
                isComplex = varargin{i};
            case "cell"
                supp = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end

    % Create the systemstruct:
    eqs = cell(s, 1);
    if ~isComplex % for real coefficients.
        for k = 1:s
            eqs{k} = [randn(size(supp{k}, 1), 1) supp{k}];
        end
    else % for complex coefficients.
        for k = 1:s
            eqs{k} = [randn(size(supp{k}, 1), 1) ...
                + 1i*randn(size(supp{k}, 1), 1) supp{k}];
        end
    end
    system = systemstruct(eqs, basis);
end