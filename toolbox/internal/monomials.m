function K = monomials(d, m, varargin)
    %MONOMIALS - Matrix with all the monomial vectors.
    %   K = MONOMIALS(d, m) creates a matrix with all the monomial vectors
    %   of the monomial basis using the graded reverse lexicographic
    %   monomial order.
    %
    %   K = MONOMIALS(..., blocksize) repeats every monomial an amount of
    %   time equal to the blocksize.
    %
    %   K = MONOMIALS(..., order) uses a user-specified monomial order.
    %
    %   Input arguments:
    %       - d (int): maximum total degree.
    %       - m (int): number of variables.
    %       - blocksize (int = 1 - optional): repetitions of each monomial.
    %       - order (function_handle = @grevlex - optional): monomial 
    %           order.
    %
    %   Output arguments:
    %       - K (int): matrix with all monomial exponent vectors.
    %
    %   See also NBMONOMIALS.
        
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.
    
    if d == 0
        K = zeros(1, m);
    else
        I = fliplr(eye(m, m));
        K = zeros(nbmonomials(d, m), m);
        K(2:m+1, :) = I;
        blockLength = ones(m, 1);
        endPoint = 1;
        for di = 2:d
            blockLength = cumsum(blockLength, "reverse");
            blockPos = [0; cumsum(blockLength)];
            endPoint = endPoint + blockLength(1);
            for ni = m:-1:1
                K(endPoint+blockPos(ni)+1:endPoint+blockPos(ni+1), :) = ...
                    K(endPoint-blockLength(ni)+1:endPoint, :) + I(ni, :);             
            end
        end
    end

    % Order the monomials according to a specific monomial ordering.
    if nargin > 2
        blocksize = 1;
        for i = 1:nargin-2
            switch class(varargin{i})
                case "function_handle"
                    order = varargin{i};
                    K = K(order(K), :);
                case "double"
                    blocksize = varargin{i};
                otherwise
                    error("Argument type is not recognized.")
            end
        end
        K = kron(K, ones(blocksize, 1));
    end
end
