function b = bkk(system)
    %BKK - Solution bound for a multivariate polynomial system.
    %   b = BKK(system) computes the Bernstein–Khovanskii–Kushnirenko (BKK)
    %   solution bound for a multivariate polynomial system.
    %
    %   Input arguments:
    %       - system (systemstruct): polynomial system.
    %
    %   Output arguments:
    %       - b (int): BKK solution bound.
    %
    %   See also BEZOUT, KUSHNIRENKO.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.
    %
    % Note: this is a naive implementation of [1, Chapter 7, Theorem 4.12].
    %
    % [1] Cox, D. A., Little, J. B., and O'Shea, D. (2005) Using Algebraic
    % Geometry, vol. 185, Springer, New York, NY, USA.

    % Verify input requirements:
    if ~isa(system, "systemstruct")
        error("Only systemstruct is allowed as input.");
    end

    % Compute the BKK solution bound via Theorem 4.12:
    b = 0;
    for k = 1:system.m
        combinations = nchoosek(1:system.m, k);
        volume = 0;
        for l = 1:size(combinations, 1)
            combo = combinations(l, :);
            B = cell(length(combo));
            for i = 1:length(combo)
                B{i} = system.supp{combo(i)};
            end
            [~, subvolume] = convhulln(minkowski(B));
            volume = volume + subvolume;
        end
        b = b + (-1)^(system.m-k)*volume;   
    end
    b = round(b);
end

% Define the required subfunction:
function result = minkowski(A)
    %MINKOWSKI - Minkowski sum of vertices.
    %   result = MINKOWSKI(A) computes the Minkowski sum of the vertices in A.
    %
    %   Input arguments:
    %       - A (cell): cell array of vertex matrices.
    %
    %   Output arguments:
    %       - result (double): vertex matrix of the Minkowski sum.

    result = A{1};
    for k = 2:length(A)
        n1 = size(result, 1);
        n2 = size(A{k}, 1);
        temp = repmat(result, n2, 1) + repelem(A{k}, n1, 1);
        result = unique(temp, "rows");
    end
end
