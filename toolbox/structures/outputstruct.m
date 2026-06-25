classdef outputstruct
    %OUTPUTSTRUCT - Class for gathering all the output details.
    %   obj = OUTPUTSTRUCT(time, nullity, shiftvalues) creates a structure 
    %   that contains the output information.
    %
    %   List of properties:
    %       - time (double): necessary computation time.
    %       - nullity (int): nullity per degree.
    %       - shiftvalues (double): numerical values of the shift polynomial.
    %       - multiplicity (logical): flag to indicate multiplicities.
    %
    %   List of methods:
    %       - disp: displays the output details.

    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    properties
        time 
        nullity
        shiftvalues
        multiplicity
    end

    methods
        function obj = outputstruct(time, nullity, shiftvalues)
            %OUTPUTSTRUCT - Constructor for the outputstruct class.
            %   obj = OUTPUTSTRUCT(time, nullity, shiftvalues) creates a
            %   structure that contains the output information.
            %
            %   Input arguments:
            %       - time (double): computation time per subalgorithm.
            %       - nullity (int): nullity per degree.
            %       - shiftvalues (double): values of the shift polynomial.
            %
            %   Output arguments:
            %       - obj (outputstruct): output structure.
            
            % Store the different output properties:
            obj.time = time;
            obj.nullity = nullity;
            obj.shiftvalues = shiftvalues;
        end 

        function disp(~)
            %DISP - Information about the output.
            %   DISP(output) displays the output details.

            fprintf("     output structure \n")
        end
    end
end
