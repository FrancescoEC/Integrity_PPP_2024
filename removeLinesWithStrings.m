function removeLinesWithStrings(inputFile, outputFile, stringList)
    % Open the input file for reading
    fidInput = fopen(inputFile, 'r');
    if fidInput == -1
        error('Cannot open input file: %s', inputFile);
    end
    
    % Open the output file for writing
    fidOutput = fopen(outputFile, 'w');
    if fidOutput == -1
        error('Cannot open output file: %s', outputFile);
    end
    
    try
        % Read and process the file line by line
        while ~feof(fidInput)
            line = fgets(fidInput); % Read a line from the input file
            if isempty(line)
                continue; % Skip empty lines
            end
            
            % Check if the first three characters match any of the strings in the list
            lineStart = strtrim(line(1:min(3, end))); % Get the first three characters
            if ~any(strcmp(lineStart, stringList))
                % If no match, write the line to the output file
                fprintf(fidOutput, '%s', line);
            end
        end
    catch ME
        % If an error occurs, close the files and rethrow the error
        fclose(fidInput);
        fclose(fidOutput);
        rethrow(ME);
    end
    
    % Close the files
    fclose(fidInput);
    fclose(fidOutput);
end