function saveFigAsPDF(figName,fontSize,overwrite)
% saveFigAsPDF Save figure as PDF file
%
% saveFigAsPDF(figName) saves the current figure in the current directory
% with the name figName, or in the path provided in figName.
%
% saveFigAsPDF(figName,fontSize) allows to set the font size (default is
% 12).
%
% saveFigAsPDF(figName,fontSize,overwrite) allows to automatically
% overwrite a saved figure with the same name when overwrite is set to TRUE.
% By default, overwrite is set to FALSE, and the program will ask whether
% you want to overwrite a saved figure with the same name.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

% set defaults if not provided
if nargin==1
    fontSize = 12;
    overwrite = false;
elseif nargin==2
    overwrite = false;
end

% set fontsize
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)

% set figure size
set(gcf,'Units','centimeters')
figPos = get(gcf,'Position');
set(gcf,'PaperPositionMode','auto','PaperUnits','centimeters',...
    'PaperSize',[figPos(3) figPos(4)]);

% check if figure exists
if ~overwrite
    if exist([figName '.pdf'],'file') == 2
        iSlash = strfind(figName,'/');
        if isempty(iSlash); iSlash = 0; end
        disp(['A figure named ' figName(iSlash(end)+1:end) ' already exists.'])
        overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
    else
        overwrite = 1;
    end
end

% save
if overwrite
    print(figName,'-dpdf');
else
    disp('Figure was not saved.')
end

end
