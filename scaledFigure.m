function h = scaledFigure(widthFactor,heightFactor)
% scaledFigure Create scaled figure window
%
% scaledFigure(widthFactor,heightFactor) creates a new figure window and
% scales its width with widthFactor and its height with heighFactor.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

h = figure;
figPos = get(h,'Position');

% scale figure size
set(h,'Position',[figPos(1:2) figPos(3)*widthFactor figPos(4)*heightFactor]);

% if figure is larger than screen, fit to screen
figPos = get(h,'Position');
screenPos = get(0,'ScreenSize');
set(h,'Position',[figPos(1:2) min([figPos(3) screenPos(3)])...
    min([figPos(4) screenPos(4)]) ]);

end
