classdef mepstruct < problemstruct
    %MEPSTRUCT   Class for multiparameter eigenvalue problems.
    %   obj = MEPSTRUCT(mat,supp) creates a problem structure based on the
    %   given coefficient matrices and support. Every coefficient matrix in 
    %   mat (= s-by-1 cell) corresponds with the monomial at the same 
    %   position in the support (= s-by-n matrix). This class inherits from 
    %   problemstruct.
    %
    %   obj = MEPSTRUCT(...,basis) also specifies the monomial basis
    %   of the given multiparameter eigenvalue problem.
    %
    %   MEPSTRUCT properties:
    %       mat - coefficient matrices.
    %       ... - properties defined in problemstruct.
    %
    %   MEPSTRUCT methods:
    %       ... - methods defined in problemstruct.  
    %
    %   See also PROBLEMSTRUCT, SYSTEMSTRUCT.

    % Copyright (c) 2024 - Christof Vermeersch

    properties
        mat
    end

    methods
        function obj = mepstruct(mat,supp,basis)
            %MEPSTRUCT   Constructor for the mepstruct class.
            %   obj = MEPSTRUCT(mat,supp) creates a problemstruct based on
            %   the given coefficient matrices and support.
            %
            %   obj = MEPSTRUCT(...,basis) also specifies the monomial 
            %   basis of the multiparameter eigenvalue problem.

            % Process the optional argument:
            if nargin < 3
                basis = "unknown";
            end

            % Create a problemstruct:
            s = numel(mat);
            [nrows,ncols] = size(mat{1});
            coef = zeros(s,nrows,ncols);
            for k = 1:s
                coef(k,:,:) = mat{k};
            end
            obj@problemstruct({coef},{supp},basis)
            obj.mat = mat;
        end
    end
end
