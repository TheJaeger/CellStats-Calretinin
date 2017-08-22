function line = getLine(imageArray)
%   Author: Siddhartha Dhiman
%   e-mail: sdhiman@buffalo.edu
%   -----------------------------------------------------------------------
%   findLine.m finds reference lines in Calretinin images and computes its
%   length and orientation in the structure array 'line'. Computation is
%   performed using Hough transform
%   -----------------------------------------------------------------------
%   Input Arguments
%       imageArray: Image file annotated with reference line
%   -----------------------------------------------------------------------
%   Output:
%       line: structural matrix containing coordinates and orientation
%             of line
%   -----------------------------------------------------------------------

%%  Extract blue channel
IB = imbinarize(imageArray(:,:,3));

%% Hough Transform
[H,T,R] = hough(IB);
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
line = houghlines(IB,T,R,P,'FillGap',5);

clearvars variables -except line
