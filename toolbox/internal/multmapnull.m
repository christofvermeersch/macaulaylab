function [A, B] = multmapnull(Z, K, G, idx, varargin)
    %MULTMAPNULL - Generalized multiplication matrices (right null space).
    %   [A, B] = MULTMAPNULL(Z, K, G, idx) sets-up the generalized
    %   multiplication matrices for a shift polynomial.
    %
    %   [A, B] = MULTMAPNULL(Z, K, G, idx) can also work with multiple 
    %   shift polynomials combined in a cell array. The resulting matrices 
    %   are then also combined in cell arrays.
    %
    %   [A, B] = MULTMAPNULL(..., blocksize) specifies the size of a block.
    %
    %   [A, B] = MULTMAPNULL(..., basis) uses a user-specified monomial
    %   basis.
    %
    %   [A, B] = MULTMAPNULL(..., order) supplements the function with
    %   a function handle for the order of the monomials.
    %
    %   Input arguments:
    %       - Z (double): (selection of) basis matrix of the null space.
    %       - K (int): monomial exponent matrix for each row of Z.
    %       - G (double): shift polynomial(s) in coefficient-exponent 
    %           format (or cell array for multiple shift polynomials).
    %       - idx (int): row indices to shift.
    %       - blocksize (int = 1 - optional): size of a block row in Z.
    %       - basis (function_handle = @monomial - optional): monomial 
    %           basis.
    %       - order (function_handle = @grevlex - optional): monomial 
    %           order.
    %
    %   Output arguments:
    %       - A (cell): first generalized multiplication matrices.
    %       - B (cell): second generalized multiplication matrices.
    %
    %   See also MULTMAPCOLUMN.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.
    
    % Process the optional parameters:
    blocksize = 1;
    basis = @monomial;
    order = @grevlex;
    for i = 1:nargin-4
        switch class(varargin{i})
            case "function_handle"
                if nargin(varargin{i}) < 2
                    order = varargin{i}; % order has one input.
                else
                    basis = varargin{i}; % basis has two inputs.
                end
            case "double"
                blocksize = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end

    % Set the correct number of shift problems:
    if iscell(G)
        isCell = true;
    else
        isCell = false;
        G = {G};
    end
    nshifts = length(G);
    m = size(G{1}, 2) - 1;

    % Compute the degree:
    d = max(sum(K, 2));

    % Determine the selection of rows in the basis:
    if blocksize == 1
        selection = order(K);
    else
        selection = (order(K(:, 1:end-1))-1)*blocksize + K(:, end);
    end

    % Define the generalized multiplication matrices:
    A = cell(nshifts, 1);
    B = cell(nshifts, 1);
    for i = 1:nshifts
        % Create first generalized multiplication map:
        A{i} = Z(idx, :);

        % Create second generalized multiplication map:
        S = rowshift(d, m, idx, G{i}, blocksize, basis, order);
        B{i} = S(1:length(idx), selection)*Z;
    end
   
    % Change output to array when input shift polynomial is an array:
    if isCell == false
        A = A{1};
        B = B{1};
    end
end
