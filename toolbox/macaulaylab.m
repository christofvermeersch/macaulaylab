function [solutions, details] = macaulaylab(problem, dend, options)
    %MACAULAYLAB - Numerical solver that uses the Macaulay matrix.
    %   solutions = MACAULAYLAB(problem) solves the problem with default
    %   options. The input and output are given as a problemstruct and 
    %   solutionstruct, respectively.
    %
    %   solutions = MACAULAYLAB(..., dend) stops the iterative procedure at
    %   the specified degree.
    %   
    %   solutions = MACAULAYLAB(..., options) uses user-specified options.
    %   The options are given in key-value pairs. For example, solutions =
    %   MACAULAYLAB(..., verbose = true).
    %
    %   solutions = MACAULAYLAB(..., dend, options) combines both the  
    %   degree and user-specified options, but since dend is a positional
    %   argument, it must be added before the option list.
    %
    %   [..., details] = MACAULAYLAB(...) also returns additional details
    %   about the algorithm and solutions.
    %
    %   Input arguments:
    %       - problem (problemstruct or subclass): polynomial problem.
    %       - dend (int = 100 - optional): last degree for the solver.
    %       - options (string - optional): key-value pairs for algorithms.
    %
    %   Output arguments:
    %       - solutions (solutionstruct): solutions and properties.
    %       - details (outputstruct): algorithm details.
    %
    %   List of options:
    %       - tol (double = 1e-10): tolerance for all the rank decisions.
    %       - clustertol (double = 1e-3): tolerance for clustering step.
    %       - polynomial (double = [randn(problem.m, 1) eye(problem.m)]):
    %           definition of the polynomial shift.
    %       - basis (function_handle = @monomial): function handle with 
    %           monomial basis.
    %       - order (function_handle = @grevlex): function handle with 
    %           monomial order.
    %       - subspace (string = "null space"): selected subspace.
    %       - blocks (string = "block-wise"): approach to deal with blocks.
    %       - enlarge (string = "sparse"): approach to enlarge subspace.
    %       - check (string = "iterative"): approach to check ranks.
    %       - jointgepsolver (string = "schur"): approach to numerically 
    %           solve joint GEP.
    %       - clustering (logical = true): flag for applying the clustering 
    %           step.
    %       - posdim (logical = false): flag for dealing with postive-
    %           dimensional solution sets.
    %       - verbose (logical = false): flag for enabling verbose output.
    %
    %   List of details:
    %       - shiftvalues (double): values of the shift polynomial.
    %       - nullity (int): nullity in the different degrees.
    %       - time (double): time required to run the algorithm.
    %       - multiplicity (logical): detection flag for multiplicities.
    %
    %   See also PROBLEMSTRUCT, SOLUTIONSTRUCT, OUTPUTSTRUCT.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.
    %   - 2026 by CV: added random Rayleigh approach for solving joint GEP.

    % Check the input and set the default options:
    arguments
        problem problemstruct
        dend double = 100
        options.tol double = 1e-10;
        options.clustertol double = 1e-3;
        options.polynomial double = [randn(problem.m, 1) eye(problem.m)];
        options.basis function_handle = @monomial;
        options.order function_handle = @grevlex;
        options.subspace string = "null space";
        options.blocks string = "block-wise";
        options.enlarge string = "sparse";
        options.check string = "iterative";
        options.jointgepsolver string = "schur";
        options.clustering logical = true;
        options.posdim logical = false;
        options.verbose logical = false;
    end

    % Parse text input:
    switch options.subspace
        case "null space"
            isNull = true;
        case "column space"
            isNull = false;
            if ~strcmp(options.enlarge, "iterative")
                options.enlarge = "recursive";
                warning("Subspace is set to enlarge recursively.")
            end
        otherwise
            error("Selected subspace invalid.")
    end
    switch options.blocks
        case "block-wise"
            isBlocks = true;
        case "row-wise"
            isBlocks = false;
        otherwise
            error("Blocks can only be block-wise or row-wise.")
    end
    switch options.check
        case "iterative"
            isRecursive = false;
        case "recursive"
            isRecursive = true;
        otherwise
            error("Check can only be iterative or recursive.")
    end

    % Fix the important variables:
    [~, dmax, m] = properties(problem);
    blocksize = size(problem.coef{1}, 3);

    % Display the properties of the problem (verbose only):
    if options.verbose
        print("about")
        print("problem", problem)
    end

    % Determine the degree of shift polynomial and start degree:
    shifts = cell(m+1, 1);
    shifts{1} = options.polynomial;
    for i = 1:m
        shifts{i+1} = [1 zeros(1, i-1) 1 zeros(1, m-i)];
    end
    dpolynomial = max(sum(options.polynomial(:, 2:end), 2));
    dstart = max(dmax, dpolynomial);
    
    % Declare necessary variables and set stabilization flag:
    time = zeros(6, 1);
    nullity = zeros(dend, 1);
    gapsize = NaN;

    % Display the properties of the algorithm (verbose only):
    if options.verbose
        print("algorithm", options.enlarge, options.subspace, options.blocks);
    end

    % Build the initial subspace:
    tic
    M = full(macaulay(problem, dstart, options.basis, options.order));
    if isNull
        Z = null(M);
        nullity(dstart) = size(Z, 2);
    else
        Z = NaN;
        r = rank(M);
        nullity(dstart) = size(M, 2) - r;
    end
    time(1) = toc;

    % Display the properties of the algorithm (verbose only):
    if options.verbose
        print("iteration", dstart, nullity(dstart), 0, gapsize)
    end

    if dstart >= dend
        error("dend (%d) must be greater than dstart (%d).", dend, dstart);
    end

    for d = dstart+1:dend
        % Enlarge the subspace:
        tic
        [M, Z] = enlarge(M, Z, problem, d, isNull, options);
        time(1) = time(1) + toc;

        % Check the subspace:
        tic
        if isNull
            nullity(d) = size(Z, 2);
        else
            r = rank(full(M));
            nullity(d) = size(M, 2) - r;
        end
        time(2) = time(2) + toc;

        if options.posdim || nullity(d) == nullity(d-dpolynomial)
            % Look for a gap zone:
            tic 
            [dgap, ma, gapsize, idx] = check(M, Z, d, m, blocksize, ...
                isNull, isBlocks, isRecursive, options);
            time(2) = time(2) + toc;

            % Break the iterative procedure when there is a gap zone:
            if gapsize >= dpolynomial
                if options.posdim
                    mb = Inf;
                else
                    mb = nullity(d);
                end
                if options.verbose
                    print("iteration", d, nullity(d), nullity(d-1), gapsize)
                end
                break
            end
        end

        % Display the properties of the algorithm (verbose only):
        if options.verbose
            print("iteration", d, nullity(d), nullity(d-1), gapsize)
        end

        % Throw error when no gap zone is found:
        if d == dend
            error("Gap zone is not found at the final degree.")
        end
    end

    % Perform a column compression (only for null space based approach):
    if isNull 
        tic
        Z = columncompr(Z, dgap + dpolynomial - 1, m, blocksize);
        time(3) = toc;
    end

    % Solve the problem:
    tic
    if strcmp(options.subspace, "null space")
        dk = dgap + dpolynomial-1;
    else
        dk = d;
    end
    K = monomials(dk, m, blocksize);
    if blocksize > 1
        K = [K kron(ones(1, nbmonomials(dk, m)), 1:blocksize)'];
    end
    switch options.subspace
        case "null space"
            if isnan(idx) % using a block-wise approach with all indices.
                idx = 1:blocksize*nbmonomials(dgap-1, m);
            end
            [A, B] = multmapnull(Z, K, shifts, idx, blocksize, ...
                options.basis, options.order);
        case "column space"
            [A, B] = multmapcolumn(M, K, shifts, idx, blocksize, ...
                options.basis, options.order);
        otherwise
            error("Solution subspace is not recognized.")
    end

    [values, shiftvalues] = jointgep(A, B, options);
    time(4) = toc;

    % Perform the optional clustering.
    if options.clustering && ma > 1
        tic
        [clusteredvalues, labels] = multiplicity([shiftvalues values], ...
            options.clustertol);
        shiftvalues = clusteredvalues(:, 1);
        values = clusteredvalues(:, 2:end);
        if max(labels) < ma
            isMultiplicity = true;
        else
            isMultiplicity = false;
        end
        time(5) = toc;
    end
 
    tic
    [maxres, res] = residuals(values, problem, options.basis);
    time(6) = toc;

    % Store the solutions and output information:
    details = outputstruct(time, nullity, shiftvalues);
    if options.clustering && ma > 1
        details.multiplicity = isMultiplicity;
    end
    solutions = solutionstruct(values);
    solutions.maxres = maxres;
    solutions.res = res;
    solutions.ma = ma;
    solutions.mb = mb;

    % Display the properties of the solution candidates (verbose only):
    if options.verbose
        print("solutions", options.clustering, solutions, details)
    end
end

% Define all the required subfunctions:
function [M, Z] = enlarge(M, Z, problem, d, isNull, options)
    %ENLARGE - Next degree of the subspaces.
    %   [M, Z] = ENLARGE(M, Z, problem, d, isNull, options) results in the
    %   next degree of the column and null space.
    %
    %   Input arguments:
    %       - M (double): current Macaulay matrix.
    %       - Z (double): current null space basis matrix.
    %       - problem (problemstruct or subclass): problem.
    %       - d (int): new degree.
    %       - isNull (logical): flag to select null space (true) or column
    %           space.
    %       - options (struct): algorithm options structure.
    %
    %   Output arguments:
    %       - M (double): updated Macaulay matrix.
    %       - Z (double): updated null space basis matrix.

    if isNull % enlarge the null space.
        switch options.enlarge
            case "iterative"
                M = full(macaulay(problem, d, options.basis, options.order));
                Z = null(M);
            case "recursive"
                rows = size(M, 1);
                M = macaulayupdate(M, problem, d, options.basis, ...
                    options.order);
                Z = nullrecrmacaulay(Z, full(M(rows+1:end, :)), options.tol);
            case "sparse"
                M = NaN;
                Z = nullsparsemacaulay(Z, problem, d, options.tol, ...
                    options.basis, options.order);
            otherwise
                error("Selected strategy to enlarge is not supported.")
        end
    else % enlarge the column space.
        switch options.enlarge
            case "iterative"
                M = full(macaulay(problem, d, options.basis, options.order));
            case "recursive"
                M = macaulayupdate(M, problem, d, options.basis, ...
                    options.order);
            otherwise
                error("Selected strategy to enlarge is not supported.")
        end
        Z = NaN;
    end
end

function [dgap, ma, gapsize, idx] = check(M, Z, d, m, vectorsize, ...
    isNull, isBlocks, isRecursive, options)
    %CHECK - Rank structure of the subspaces.
    %   [...] = CHECK(...) checks the rank structure of the subspaces.
    %
    %   Input arguments:
    %       - M (double): Macaulay matrix.
    %       - Z (double): null space basis matrix.
    %       - d (int): current degree.
    %       - m (int): number of variables.
    %       - vectorsize (int): block size (eigenvector length).
    %       - isNull (logical): flag for null space (true) or column space.
    %       - isBlocks (logical): flag for block-wise (true) or row-wise 
    %           approach.
    %       - isRecursive (logical): flag for recursive (true) or iterative
    %           check.
    %       - options (struct): algorithm options structure.
    %
    %   Output arguments:
    %       - dgap (int): first degree of the gap zone (NaN if no gap zone).
    %       - ma (int): number of affine solutions.
    %       - gapsize (int): size of the gap zone (in number of degree blocks).
    %       - idx (int): indices of the standard monomials.

    if isBlocks && isNull % find gap via blocked or non-blocked approach.
            % Set the standard monomials to NaN:
            idx = NaN;
            idxinf = NaN;
            [dgap, ma, gapsize] = gap(Z, m, vectorsize, isRecursive, ...
                options.tol);
    else
        if isNull
            % Find the standard monomials:
            idx = stdmonomials(Z, options.tol);
        else
            % Find the standard monomials:
            sinv = stdmonomials(full(flip(M')), options.tol);
            idx = setdiff(1:size(M, 2), size(M, 2) - sinv + 1);
        end
        % Determine the number of affine solutions:
        [dgap, ma, gapsize] = gapstdmonomials(idx, d, m, vectorsize);
        idx = idx(1:ma);
    end
end

function print(component, varargin)
    %PRINT - Output of MacaulayLab via fprintf.
    %   PRINT(component) generates the output of MacaulayLab for a certain
    %   component of the solver. It is an auxiliary function that uses
    %   the built-in fprintf.
    %
    %   PRINT(..., varargin) also processes the additional information in
    %   the variable argument list.
    %
    %   Input arguments:
    %       - component (string): string identifying the output component
    %           "about" | "problem" | "algorithm" | "iteration" | 
    %           "solutions".
    %       - varargin (cell): additional data for the selected component.

    LINEWIDTH = 75;
    linestr = pad("", LINEWIDTH, "_");

    switch component
        case "about"
            fileId = fopen("about.txt");
            about = fscanf(fileId, "%c");
            fclose(fileId);
            % fprintf(linestr);
            fprintf("\n");
            fprintf(about);
            fprintf("\n");
            fprintf(linestr);
            fprintf("\n");

        case "problem"
            fprintf("\n");
            fprintf(pad("PROBLEM INFO", LINEWIDTH, "both"));
            fprintf("\n");
            info(varargin{1});
            fprintf(linestr);
            fprintf("\n"); 

        case "algorithm"
            fprintf("\n");
            fprintf(pad("ALGORITHM OUTPUT", LINEWIDTH, "both"));
            fprintf("\n");
            fprintf("    using a Macaulay matrix algorithm")
            fprintf(" with a %s construction of \n", varargin{1});
            fprintf("    the %s and %s shifts \n \n", varargin{2}, ...
                varargin{3});
            fprintf("    | degree        | nullity        | increase" ...
                + "       | gap zone      | \n")
            fprintf("    |------------------------------------------" + ...
            "-----------------------| \n")
          
        case "iteration"
            degree = sprintf("%d", varargin{1});
            degree = degree + pad("", 13-strlength(degree), " ");
            nullity = sprintf("%d", varargin{2});
            nullity = nullity + pad("", 14-strlength(nullity), " ");
            increase = sprintf("%d", varargin{2}-varargin{3});
            increase = increase + pad("", 14-strlength(increase), " ");
            gapzone = sprintf("%d", varargin{4});
            gapzone = gapzone + pad("", 13-strlength(gapzone), " ");
            fprintf("    | %s | %s | %s | %s | \n", ...
            degree, nullity, increase, gapzone)

        case "solutions"
            fprintf(linestr);
             fprintf("\n");
            fprintf("\n");
            fprintf(pad("SOLUTIONS", LINEWIDTH, "both"));
            fprintf("\n");
            info(varargin{2});
            if varargin{1}
                if varargin{3}.multiplicity
                    fprintf("    multiple solutions detected \n");
                else
                    fprintf("    no multiple solutions detected \n");
                end
            end
            fprintf("\n")
            fprintf("    total computation is %.1e seconds \n", ...
                sum(varargin{3}.time))
            fprintf("    maximum residual is %.1e \n", varargin{2}.maxres)
            fprintf(linestr);
            fprintf("\n"); 
    end 
end

function [values, shiftvalues] = jointgep(A, B, options)
    %JOINTGEP - Compute the eigenvalues of a joint GEP.
    %   [values, shiftvalues] = JOINTGEP(A, B, options) computes the
    %   eigenvalues of the joint GEP determined in the cells of A and B.
    %   Two approaches, i.e., the Schur factorization and random Rayleigh
    %   technique, are supported.
    %
    %   Input arguments:
    %       - A (cell array): first generalized multiplication matrices.
    %       - B (cell array): second generalized multiplication matrices.
    %       - options: algorithm options (constructed by main function).
    %
    %   Output arguments:
    %       - values (double): matrix of solution values.
    %       - shiftvalues (double): values of the shift polynomial.

    if strcmp(options.jointgepsolver, "schur")
        [Q, S] = schur(pinv(full(A{1}))*full(B{1}), "complex");
        shiftvalues = diag(S);
        values = zeros(length(S), length(A) - 1);

        for i = 1:length(A) - 1
            S = Q'*pinv(full(A{i+1}))*B{i+1}*Q;
            values(:, i) = diag(S);
        end
    elseif strcmp(options.jointgepsolver, "randomrayleigh")
        if ~strcmp(options.subspace, "null space")
            error("Random Rayleigh approach only works in the null space.")
        end
        
        [X, S, Y] = eig(pinv(full(A{1}))*full(B{1}));
        shiftvalues = diag(S);
        values = zeros(length(S), length(A) - 1);
    
        den = zeros(size(A{1}, 1), 1);
        Gamma = A{1}*X;
        for i = 1:size(A{1}, 1)
            den(i) = Y(:,i)'*Gamma(:,i);
        end
        for i = 1:length(A) - 1
            Gamma = B{i + 1}*X;
            for j = 1:size(A{1}, 1)
                values(j, i) = (Y(:, j)'*Gamma(:, j))/den(j);
            end
        end
    else
        error("Solver for joint GEP does not exist.")
    end
end


