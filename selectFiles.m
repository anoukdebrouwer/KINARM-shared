function selectedFiles = selectFiles(name,filesOrFolders)
% selectFiles  Select files or folders from a list
%
% selectedFiles = selectFiles('name') prints all files and folders that
% match the string names, waits for user input to select files, and returns
% a struct with file info of the selected files. name can be a file or
% folder in the current directory, or name can specify an absolute or
% relative path. name can include wildcards (*), e.g., 's*' or '*.mat'.
%
% selectedFiles = selectFiles('name','kind') allows to specify to look for
% files or folders only. The 'kind' argument may be 'files' or 'folders'.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

if nargin==1
    filesOrFolders = 'all';
end

% find folders or files
files = dir(name);
firstLetter = cellfun(@(x) x(1),{files.name}');
files = files(isstrprop(firstLetter,'alphanum')); % discard hidden files on Mac
if strcmp(filesOrFolders,'files')
    files = files([files.isdir]==0);
elseif strcmp(filesOrFolders,'folders')
    files = files([files.isdir]==1);
end
fileNames(:,1) = num2cell((1:length(files))');
fileNames(:,2) = {files.name}' % show file numbers and names

% ask for user input
if ~isempty(fileNames)
    fileNumbers = input('Select data files or folders by typing the numbers \n(e.g., [1:5,7]) or press enter to select all: ');
    if fileNumbers>0
        selectedFiles = files(fileNumbers);
    else
        selectedFiles = files;
    end
else
    selectedFiles = [];
    disp('No folders or files found')
end

end
