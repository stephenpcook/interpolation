function out = interp1cubicL(x, M, y)
%INTERP1CUBICL Interpolation using cubic lagrange interpolation
%
% pp = INTERP1CUBICL(x, M) creates a piecewise cubic lagrange interpolant
% defined as f(x(i)) = M(i), and returns a piecewise polynomial structure
% for use with e.g. ppval. The end intervals are defined as linear, whilst
% the internal intervals are defined with the cubic lagrange interpolant
% for the 4 surrounding points.
%
% M_y = INTERP1CUBICL(x, M, y) is the same as above, but returns the
% piecewise polynomial evaluated at y.
%
% Example: If X = (0:5) and M = sin(X), then we can make and plot the pp
%
%     pp = interp1cubicL(X,M);
%     y = 0:0.1:5;
%     plot(y,ppval(pp,y),'-',X,M,'+');
%
% or directly with
%
%     plot(y, interp1cubicL(X,M,y),'-',X,M,'+');
%
% See also: INTERP3LIM INTERP_ENO, PPVAL_LIM

% Author: Stephen P. Cook <s.cook@bath.ac.uk>
% Date: 11-02-2016

x_N = length(x);
% Check about strict monotonicity of x and x_N>2 ?
% Check M same length as x?
% Check all y lie in the range of x?

% Set up the pointwise polynomial matrix!
pp_mat = zeros((x_N-1),4);
% These may all be upside down.
% Linear in first and last interval
pp_mat(1,:) = [0,0,(M(2)-M(1))/(x(2)-x(1)),M(1)];
pp_mat(x_N-1,:) = [0;0;(M(end)-M(end-1))/(x(2)-x(1));M(end-1)];
for ii=2:(x_N-2)
    xtilde = [x((ii-1):(ii+2))] - x(ii);
    m_scale = [...
        (xtilde(1)-xtilde(2))*(xtilde(1)-xtilde(3))*(xtilde(1)-xtilde(4));...
        (xtilde(2)-xtilde(1))*(xtilde(2)-xtilde(3))*(xtilde(2)-xtilde(4));...
        (xtilde(3)-xtilde(1))*(xtilde(3)-xtilde(2))*(xtilde(3)-xtilde(4));...
        (xtilde(4)-xtilde(1))*(xtilde(4)-xtilde(2))*(xtilde(4)-xtilde(3));...
        ];
    mtilde = [M((ii-1):(ii+2))];
    mtilde = mtilde(:)./m_scale;
    A = [1, 1, 1, 1;...
        -(xtilde(2)+xtilde(3)+xtilde(4)), -(xtilde(1)+xtilde(3)+xtilde(4)),...
        -(xtilde(1)+xtilde(2)+xtilde(4)), -(xtilde(1)+xtilde(2)+xtilde(3));...
        xtilde(2)*xtilde(3)+xtilde(2)*xtilde(4)+xtilde(3)*xtilde(4),...
        xtilde(1)*xtilde(3)+xtilde(1)*xtilde(4)+xtilde(3)*xtilde(4),...
        xtilde(1)*xtilde(2)+xtilde(1)*xtilde(4)+xtilde(2)*xtilde(4),...
        xtilde(1)*xtilde(2)+xtilde(1)*xtilde(3)+xtilde(2)*xtilde(3);...
        -xtilde(2)*xtilde(3)*xtilde(4), -xtilde(1)*xtilde(3)*xtilde(4),...
        -xtilde(1)*xtilde(2)*xtilde(4), -xtilde(1)*xtilde(2)*xtilde(3)];
    pp_mat(ii,:) = (A*mtilde(:))';

end % for loop
% pp format
if nargin==2
    out= mkpp(x,pp_mat);
elseif nargin==3
    out = ppval(mkpp(x,pp_mat),y);
end % function
