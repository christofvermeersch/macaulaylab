function tensor = tensorization(mat)
    %TENSORIZATION - Conversion from matrix to tensor representation.
    %   tensor = TENSORIZATION(coef) creates a tensor of coefficient
    %   matrices from the cell array with coefficient matrices.
    %
    %   Input arguments:
    %       - mat (cell): coefficient matrices.
    %
    %   Output arguments:
    %       - tensor (double): coefficient tensor.
    %
    %   See also MEPSTRUCT, PROBLEMSTRUCT.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    s = numel(mat);
    [nrows, ncols] = size(mat{1});
    tensor = permute(reshape(horzcat(mat{:}), [nrows, ncols, s]), ...
        [3, 1, 2]);
end