function [theta,phi,r] = cart2sphd(x,y,z,maxAngle)
% cart2sphd Convert xyz position to spherical coordinates in degrees
%
% [theta,phi,r] = cart2sphd(x,y,z) converts xyz position to spherical
% coordinates with angles theta (azimuth) and phi (elevation) in degrees.
% Theta = 0 degrees is 'north' and clockwise rotations are positive. For
% standard spherical coordinates, use MATLAB's built-in cart2sph function.
%
% [theta,phi,r] = cart2sphd(x,y,z,180) returns angles theta in the range of
% -180 to 180 degrees.

% MIT License
% Copyright (c) 2021 Anouk de Brouwer

% set default maximum angle angle if unspecified
if nargin==2
    maxAngle = 360;
end

% calculate distance
r = calc3Ddistance([0 0 0],[x(:) y(:) z(:)]);

% calculate azimuth
theta = xy2compass(x,y,maxAngle);

% calculate elevation
phi = acosd(z./r(:));