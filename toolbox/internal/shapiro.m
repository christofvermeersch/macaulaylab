function b = shapiro(mep)
    %SHAPIRO - Solution bound for a multiparameter eigenvalue problem.
    %   b = SHAPIRO(mep) computes the Shapiro-Shapiro solution bound
    %   for a linear multiparameter eigenvalue problem.
    %
    %   Input arguments:
    %       - mep (mepstruct): linear multiparameter eigenvalue problem.
    %
    %   Output arguments:
    %       - b (int): Shapiro-Shapiro solution bound.
    %
    %   See also HKP.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Verify input requirements:
    if ~isa(mep, "mepstruct")
        error("Only mepstruct is allowed as input.");
    end

    % Compute Shapiro-Shapiro number:
    b = nchoosek(mep.m+mep.l-1, mep.m);
end 
