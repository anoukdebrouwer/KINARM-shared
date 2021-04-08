function [Index]=getColumnIndex(ColumnHeaders,VariableName)
% getColumnIndex  Get index of data column by variable name

% Jason Gallivan
% Modified by Anouk de Brouwer

NumColumns = size(ColumnHeaders,2);

i = 1;
ColumnNumber = 1;
Done = 0;
while Done==0
    String = char(ColumnHeaders(ColumnNumber));
    if strcmpi(String,VariableName)
        % found it
        Index(i) = ColumnNumber;
        i = i+1;
    end
    if Done == 0
        ColumnNumber = ColumnNumber + 1;
        if ColumnNumber > NumColumns
            % not found
            Done = 1;
            if ~exist('Index','var')
                Index = 0;
            end
        end
    end
end