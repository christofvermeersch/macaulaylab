function output = information(input)
    %INFORMATION - Provide information about MacaulayLab.
    %   INFORMATION(input) prints the requested information. The input
    %   parameter is a string, describing the requested information.
    %   
    %   output = INFORMATION(...) returns the information as an output
    %   string.
    %
    %   Input arguments:
    %       - input (string): string identifying the requested information
    %           ("version" | "license").
    %
    %   Output arguments:
    %       - output (string): string with the requested information.
    
    % Copyright (c) 2026 - Christof Vermeersch
    %
    % Updates:
    %   - 2026 by CV: updated documentation and comments.

    switch input
        case "version"
            text = "1.1.0 \n";
        case "license"
            fileId = fopen("license.txt");
            text = fscanf(fileId, "%c");
            fclose(fileId);
        otherwise
            error("Requested information is not available.");
    end

    if nargout > 0
        output = text;
    else
        fprintf(text);
    end
end