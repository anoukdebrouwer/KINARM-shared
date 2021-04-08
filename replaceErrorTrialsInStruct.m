function S = replaceErrorTrialsInStruct(S,out,replaceZerosWithNans)
% replaceErrorTrialsInStruct  Replace the values in a data struct with NaNs
% for trials categorized as error trials
%
% S = replaceErrorTrialsInStruct(S,out,replaceZerosWithNans)

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

%% Replace zeros and error trials by NaNs

names = fieldnames(S); % field names
names(strcmp(names,'description')) = []; % remove field with description

for f = 1 : length(names) % loop over struct fields
    currField = S.(names{f});
    [nRows,nCols,nStacks] = size(currField);
    if nRows==size(out,1) % check if number of rows is identical
        
        if nCols==size(out,2) % check if number of columns is identical
            if isnumeric(currField) % numeric field
                % replace zeros with NaNs
                if replaceZerosWithNans
                    currField(currField==0) = NaN;
                end
                % replace error trials with NaNs
                for i = 1 : nStacks
                    currField_temp = currField(:,:,i);
                    currField_temp(out) = NaN;
                    currField(:,:,i) = currField_temp;
                end
                
            elseif iscell(currField) % cell field
                % replace empty cells with NaNs
                if replaceZerosWithNans
                    emptyCell = cellfun(@isempty,currField);
                    currField(emptyCell) = {NaN};
                end
                % replace error trials with NaNs
                for i = 1 : nStacks
                    currField_temp = currField(:,:,i);
                    currField_temp(out) = {NaN};
                    currField(:,:,i) = currField_temp;
                end
                
            end
            S.(names{f}) = currField;
            
        elseif size(out,2)==1 && nCols>1 % out is vector, field is matrix
            if isnumeric(currField) % numeric field
                % replace zeros with NaNs
                if replaceZerosWithNans
                    currField(currField==0) = NaN;
                end
                % replace error trials with NaNs
                currField(out,:) = NaN;
                
            elseif iscell(currField) % cell field
                % replace error trials with NaNs
                currField(out,:) = {NaN};
                
            end
            S.(names{f}) = currField;
            
        else
            keyboard
        end
        
    else
        disp(['Number of rows is not identical. No error trials removed from field ' names{f}])
    end
end
