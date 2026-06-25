classdef mepstruct < problemstruct
    %MEPSTRUCT - Class for multiparameter eigenvalue problems.
    %   obj = MEPSTRUCT(mat, supp) creates a problem structure based on the
    %   given coefficient matrices and support. This class inherits from
    %   problemstruct.
    %
    %   obj = MEPSTRUCT(..., basis) also specifies the monomial basis
    %   of the given multiparameter eigenvalue problem.
    %
    %   List of properties:
    %       - mat (cell): coefficient matrices.
    %       - ...: properties defined in problemstruct.
    %
    %   List of methods:
    %       - ...: methods defined in problemstruct.  
    %
    %   See also PROBLEMSTRUCT, SYSTEMSTRUCT.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    properties
        mat
    end

    methods
        function obj = mepstruct(mat, supp, basis)
            %MEPSTRUCT - Constructor for the mepstruct class.
            %   obj = MEPSTRUCT(mat,supp) creates a problemstruct based on
            %   the given coefficient matrices and support.
            %
            %   obj = MEPSTRUCT(...,basis) also specifies the monomial
            %   basis of the multiparameter eigenvalue problem.
            %
            %   Input arguments:
            %       - mat: cell array of coefficient matrices (s-by-1).
            %       - supp: support matrix of monomial exponents (s-by-m).
            %       - basis (function_handle - optional): monomial basis.
            %
            %   Output arguments:
            %       - obj (mepstruct): multiparameter eigenvalue problem.

            % Process the optional argument:
            if nargin < 3
                basis = NaN;
            end

            % Create a problemstruct:
            coef = tensorization(mat);
            obj@problemstruct({coef}, {supp}, basis)
            obj.mat = mat;
        end
    end
end
