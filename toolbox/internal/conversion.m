function mep = conversion(system, eigenvectorIdx)
    %CONVERSION - Multiparameter eigenvalue problem conversion.
    %   mep = CONVERSION(system, eigenvectorIdx) gives in a multiparameter
    %   eigenvalue problem based on the given system and the selection of
    %   some variables as the elements of the eigenvector.
    %
    %   Input arguments:
    %       - system (systemstruct): polynomial system.
    %       - eigenvectorIdx (int): indices of the variables that serve as
    %           eigenvector elements.
    %
    %   Output arguments:
    %       - mep (mepstruct): multiparameter eigenvalue problem.
    %
    %   See also MEPSTRUCT, SYSTEMSTRUCT, PROBLEMSTRUCT.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.
    %
    % This conversion code is an adaptation of a script written by Sibren
    % Lagauw.

    % Verify input requirements:
    if ~isa(system, "systemstruct")
        error("Only systemstruct is allowed as input.");
    end

    % Divide the variables into parameters and eigenvector elements:
    allIdx = 1:system.m;
    paramIdx = setdiff(allIdx, eigenvectorIdx);

    % Compute the degree of the multiparameter eigenvalue problem:
    dmep = 0;
    for i = 1:system.s
        M = system.eqs{i};
        dmep = max(dmep, max(sum(M(:, paramIdx+1), 2)));
    end

    % Create matrices with monomials:
    paramMonomials = monomials(dmep, length(paramIdx));
    eigenvectorMonomials = monomials(1, length(eigenvectorIdx));

    % Compose coefficient matrices: 
    mepmat = cell(nbmonomials(dmep, length(paramIdx)), 1);
    for i = 1:size(paramMonomials, 1) % loop over the different matrices.
        M = zeros(system.s, length(eigenvectorIdx)+1);
        for j = 1:system.s % loop over the different polynomials.
            poly = system.eqs{j};
            for k = 1:size(poly, 1) % loop over every coefficient.
                tmp = poly(k, paramIdx+1);
                % Check if the coefficient belongs to the matrix:
                if isequal(paramMonomials(i, :), tmp)
                    % Find correct column of matrix: 
                    monomial = poly(k, eigenvectorIdx+1);
                    [~, idx] = ismember(monomial, eigenvectorMonomials, ...
                        "rows");
                    if idx > 0
                        M(j, idx) = poly(k, 1);
                    end
                end
            end
        end
        mepmat{i} = M;
    end

    % Create the mepstruct:
    mep = mepstruct(mepmat, paramMonomials);
end