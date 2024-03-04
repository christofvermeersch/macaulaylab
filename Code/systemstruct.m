classdef systemstruct < problemstruct
    properties
        eqs
    end

    methods
        function obj = systemstruct(eqs,basis)
            if ~exist('basis','var')
                basis = 'monomial';
            end
            s = numel(eqs);
            coef = cell(s,1);
            supp = cell(s,1);
            for k = 1:s
                polynomial = eqs{k};
                coef{k} = polynomial(:,1);
                supp{k} = polynomial(:,2:end);
            end
            obj@problemstruct(coef,supp,basis)
            obj.type = 'system';
            obj.eqs = eqs;
        end 

    end
end