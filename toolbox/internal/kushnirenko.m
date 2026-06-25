function b = kushnirenko(system)
    %KUSHNIRENKO - Solution bound for a multivariate polynomial system.
    %   b = KUSHNIRENKO(system) computes the Kushnirenko bound for a
    %   multivariate polynomial system.
    %
    %   Input arguments:
    %       - system (systemstruct): polynomial system.
    %
    %   Output arguments:
    %       - b (int): Kushnirenko solution bound.
    %
    %   See also BEZOUT, BKK.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Verify input requirements:
    if ~isa(system, "systemstruct")
        error("Only systemstruct is allowed as input.");
    end

    % Determine the unique support of the system:
    uniqueSupport = unique(cat(1, system.supp{:}), "rows");

    % Compute the Kushnirenko solution bound:
    [~, volume] = convhulln(uniqueSupport);
    b = round(factorial(system.m)*volume);
end