function theta = xy2compass(x,y,maxAngle)
% xy2compass  Convert xy position to compass angle in degrees
%
% theta = xy2compass(x,y) converts xy position to angle in degrees
% where 0 degrees is north and clockwise rotations are positive.
% Default values are in the range of 0 to 360 degrees.
%
% theta = xy2compass(x,y,180) returns angles in the range of -180 to
% 180 degrees.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

% set default angle if unspecified
if nargin==2
    maxAngle = 360;
end

% compute inverse tangens in degrees
theta = atan2d(y,x);
theta = -theta +90; % unit circle to compass angles

% convert to [0 360] unless a maximum angle of 180 is specified
if maxAngle~=180
    neg = theta<0;
    theta(neg) = theta(neg)+360;
end
