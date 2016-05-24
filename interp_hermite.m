function [Vq,Dq] = interp_hermite(X,V,Xq,d_type)
%INTERP_HERMITE 1D Hermite interpolation
%
% [Vq,Dq] = INTERP_HERMITE(X,V,Xq,d_type) evaluates the Hermite interpolant
% Vq and its derivative Dq for a function V=F(X) at query points Xq using
% gradient estimate determined by the string d_type passed into
% calc_gradients (default 'hyman').
%
% [Vq,Dq] = INTERP_HERMITE(X,V,Xq) same as above with d_type='hyman'.
%
% X and V are vectors of the same length, and Vq and Dq will be returned
% with the same shape as Xq. For any points in Xq which are also in X, the
% corresponding Dq is the derivative from the right.
%
% Example: Here we show the interpolant of a sine wave with hermite
% interpolation with the default gradient estimator.
%     x = 0:6;
%     y = sin(x);
%     xx = linspace(0,6);
%     [yy,dd] = interp_hermite(x,y,xx);
%     plot(xx,yy,'b-',x,y,'r-')
%
% See also: CALC_GRADIENTS INTERP_HERMITE_LIM EVAL_HERMITE

% Uses external functions calc_gradients and eval_hermite.

if ( nargin==3 )
    d_type = 'hyman';
end % if

X = X(:);
V = V(:);

if length(X)<4
  error('Requires at least 4 data points for X and V.');
end % if
if any(Xq<min(X) | Xq>max(X))
  error('Query points in Xq lie outside the range of X.');
end % if

D = calc_gradients(X,V,d_type);

[Vq,Dq] = eval_hermite(X,V,D,Xq);

end % function interp_hermite