function distance = calc3Ddistance(xyz1,xyz2)
% calc3Ddistance  Calculate the Euclidean distance between two points
%
% distance = calc3Ddistance(xyz1) calculates the length of vector xyz1
% (from origin [0 0 0]), where xyz1 can be a single point or a set of
% points in space.
%
% distance = calc3Ddistance(xyz1,xyz2) calculates the distance between xyz1
% and xyz2, where xyz2 can be a single point or a set of points matching
% the number of points in xyz1.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

if nargin==1            % if xyz2 is undefined, set to zero
    xyz2 = [0 0 0];
end
n1 = size(xyz1,1);      % number of coordinate sets
n2 = size(xyz2,1);

distance = zeros(max([n1 n2]),1);
if n1>1 && n2==1        % multiple xyz1 sets, one xyz2 set
    for i = 1 : n1
        x2 = (xyz1(i,1) - xyz2(1,1))^2;
        y2 = (xyz1(i,2) - xyz2(1,2))^2;
        z2 = (xyz1(i,3) - xyz2(1,3))^2;
        distance(i) = sqrt(x2+y2+z2);
    end
elseif n1==1 && n2>1    % one xyz1 set, multiple xyz2 sets
    for i = 1 : n2
        x2 = (xyz1(1,1) - xyz2(i,1))^2;
        y2 = (xyz1(1,2) - xyz2(i,2))^2;
        z2 = (xyz1(1,3) - xyz2(i,3))^2;
        distance(i) = sqrt(x2+y2+z2);
    end
elseif n1==n2           % same number of sets
    for i = 1 : n1
        x2 = (xyz1(i,1) - xyz2(i,1))^2;
        y2 = (xyz1(i,2) - xyz2(i,2))^2;
        z2 = (xyz1(i,3) - xyz2(i,3))^2;
        distance(i) = sqrt(x2+y2+z2);
    end
end
