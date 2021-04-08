function [i,iStartEnd] = findIntervals(X,n)
% findIntervals  Find start and end indices of intervals during which one
% or multiple conditions are nonzero
%
% i = findIntervals(X,n) returns the first index of when all columns in X
% are nonzero for n or more rows. X should be a vector or a matrix with the
% conditions in columns.
%
% [i,iStartEnd] = findIntervals(X,n) returns a k-by-2 matrix with the start
% and end row indices of all k intervals in which all columns in X are
% nonzero for n or more rows.
%
% If X contains no nonzero elements, findIntervals returns NaN.
%
% Example use case: find when both the x and y position of the hand are in
% a target area for at least 10 measurement samples.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

% make sure the input is in columns
if size(X,1)<size(X,2)
    X = X';
end

% find start and end indices of when all conditions are fulfilled
allXTrue = prod(X,2);
startTrue = find(diff(allXTrue)==1);
endTrue = find(diff(allXTrue)==-1);

% make sure the number of start and end indices is the same
% if condition is true from first row
if allXTrue(1) && (isempty(startTrue) || (startTrue(1)>endTrue(1)))
    startTrue = [0; startTrue];
end
% if condition is true until last row
if allXTrue(end) && (isempty(endTrue) || (startTrue(end)>endTrue(end)))
    endTrue = [endTrue; length(allXTrue)];
end

% check length of intervals during which all conditions are true
trueStartEnd = [startTrue endTrue];
if ~isempty(trueStartEnd)
    nTrue = trueStartEnd(:,2) - trueStartEnd(:,1);
    trueStartEnd = trueStartEnd(nTrue>=n,:);
    if ~isempty(trueStartEnd)
        i = trueStartEnd(1)+1;
        iStartEnd = [trueStartEnd(:,1)+1 trueStartEnd(:,2)];
    else
        i = NaN;
        iStartEnd = NaN;
    end
else
    i = NaN;
    iStartEnd = NaN;
end
