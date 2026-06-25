classdef systemstruct < problemstruct
    %SYSTEMSTRUCT - Class for multivariate polynomial systems.
    %   obj = SYSTEMSTRUCT(eqs) creates a problem structure based on the
    %   given polynomial equations. The polynomials can also be symbolic 
    %   expressions. This class inherits from problemstruct.
    %
    %   obj = SYSTEMSTRUCT(..., basis) also specifies the monomial basis
    %   of the multivariate polynomial system.
    %
    %   List of properties:
    %       - eqs (cell | string): polynomial equations.
    %       - map (cell): mapping of symbolic variables to the indices.
    %       - ...: properties defined in problemstruct.
    %
    %   List of methods:
    %       - ...: methods defined in problemstruct.   
    %
    %   See also PROBLEMSTRUCT, MEPSTRUCT.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2025 by CV: added support for univariate symbolic polynomials.
    %   - 2026 by CV: updated documentation and comments.
    %
    % The included conversion from a symbolic representation to the
    % system structure used in MacaulayLab is based on the code of Oscar
    % Mauricio Agudelo.
    

    properties
        eqs
        map
    end

    methods
        function obj = systemstruct(eqs, basis)
            %SYSTEMSTRUCT - Constructor for the systemstruct class.
            %   obj = SYSTEMSTRUCT(eqs) creates a problemstruct based on
            %   the given multivariate polynomial equations. The
            %   polynomials can also be symbolic expressions.
            %
            %   obj = SYSTEMSTRUCT(..., basis) also specifies the monomial
            %   basis of the multivariate polynomial system.
            %
            %   Input arguments:
            %       - eqs (cell | string): polynomial equations.
            %       - basis (function_handle - optional): monomial basis.
            %
            %   Output arguments:
            %       - obj (systemstruct): polynomial system structure.

            % Process the optional argument:
            if nargin < 2
                basis = NaN;
            end

            % Convert symbolic expressions into numerical polyonomials:
            s = numel(eqs);
            variables = NaN;
            if isa(eqs, "sym") || isa(eqs{1}, "sym")
                if isa(eqs, "cell")
                    eqs = cell2sym(eqs);
                end
                    
                symboliceqs = eqs;
                variables = symvar(symboliceqs);
                if isscalar(variables)
                    syms additionalvariable 
                    symboliceqs(s+1) = additionalvariable;
                    variables = symvar(symboliceqs);
                    idx = find(variables == additionalvariable);
                end
                variablenames = char(variables);
                eqs = cell(s, 1);
                for i = 1:s
                    t = feval(symengine, "poly2list", ...
                        vpa(symboliceqs(i)), variablenames);
                    poly = zeros(length(t), length(variables)+1);
                    for j = 1:length(t)
                        monomial = t(j);
                        poly(j, 1) = monomial(1);
                        poly(j, 2:end) = monomial(2);
                    end
                    if exist("additionalvariable", "var")
                        eqs{i} = [poly(:, 1:idx) poly(:, idx+2:end)];
                    else
                        eqs{i} = poly;
                    end
                end
            end
            
            % Create a problemstruct:
            coef = cell(s, 1);
            supp = cell(s, 1);
            for k = 1:s
                polynomial = eqs{k};
                coef{k} = polynomial(:, 1);
                supp{k} = polynomial(:, 2:end);
            end
            obj@problemstruct(coef, supp, basis)
            obj.eqs = eqs;
            obj.map = variables;
        end 
    end
end
