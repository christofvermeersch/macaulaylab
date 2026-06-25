function S = rowshift(d, m, idx, poly, varargin)
    %ROWSHIFT - Matrix of shifts.
    %   S = ROWSHIFT(d, m, idx, poly) creates a sparse shift matrix in 
    %   which rows with given indices are shifted by a certain shift 
    %   polynomial.
    %
    %   S = ROWSHIFT(..., blocksize) considers problems with eigenvectors.
    %
    %   S = ROWSHIFT(..., basis) performs the shift in a user-specified
    %   monomial basis.
    %
    %   S = ROWSHIFT(..., order) performs the shift in a user-specified
    %   monomial order.
    %
    %   Input arguments:
    %       - d (int): degree of the Macaulay matrix.
    %       - m (int): number of variables.
    %       - idx (int): row indices to shift.
    %       - poly (double): shift polynomial in coefficient-exponent 
    %           format.
    %       - blocksize (int = 1 - optional): eigenvector length.
    %       - basis (function_handle = @monomial - optional): monomial 
    %           basis.
    %       - order (function_handle = @grevlex - optional): monomial 
    %           order.
    %
    %   Output arguments:
    %       - S (double): sparse shift matrix.
    %
    %   See also ROWRECOMB.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
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

    % Determine basis dilation:
    bd = size(basis(ones(1, m), ones(1, m)), 1);

    % Create extended monomials matrix:
    K = [monomials(d, m, order, blocksize) kron(ones(1, ...
        nbmonomials(d, m)), 1:blocksize)'];

    % Determine shifts of all the rows with given indices:
    shifts = basis(K(idx, 1:end-1), poly(:, 2:end));
    positions = order(shifts(:, 2:end));

    % Construct the matrix of shifts:
    rowIdx = kron(1:length(idx), ones(1, size(poly, 1)*bd)).';
    colIdx = (positions-1)*blocksize+kron(K(idx, end), ...
        ones(size(poly, 1)*bd, 1));
    values = shifts(:, 1).*kron(ones(bd*length(idx), 1), poly(:, 1));
    S = sparse(rowIdx, colIdx, values, length(idx), ...
        nbmonomials(d, m, blocksize));
end
