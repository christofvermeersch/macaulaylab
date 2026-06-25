function b = hkp(mep)
    %HKP - Solution bound for a multiparameter eigenvalue problem.
    %   b = HKP(mep) computes the Hochstenbach-Kosir-Plestenjak (HKP)
    %   solution bound for a polynomial multiparameter eigenvalue problem.
    %
    %   Input arguments:
    %       - mep (mepstruct): multiparameter eigenvalue problem.
    %
    %   Output arguments:
    %       - b (int): HKP solution bound.
    %
    %   See also SHAPIRO.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    % Verify input requirements:
    if ~isa(mep, "mepstruct")
        error("Only mepstruct is allowed as input.");
    end
    
    % Compute the HKP number:
    b = mep.dmax^mep.m*nchoosek(mep.m+mep.l-1, mep.m);
end
