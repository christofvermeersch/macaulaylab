function b = bezout(system)
    %BEZOUT - Bézout number for a multivariate polynomial system.
    %   b = BEZOUT(system) computes the Bezout number for a
    %   multivariate polynomial system.
    %
    %   Input arguments:
    %       - system (systemstruct or problemstruct): polynomial system.
    %
    %   Output arguments:
    %       - b (int): Bézout number.
    %
    %   See also KUSHNIRENKO, BKK.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Verify input requirements:
    if ~isa(system, "systemstruct") && ~isa(system, "problemstruct")
        error("Only systemstruct is allowed as input.");
    end
    
    % Compute Bézout number.
    b = prod(vertcat(system.di{:}));
end