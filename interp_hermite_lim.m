function [Vq,Dq] = interp_hermite_lim(X,V,Xq,d_type)
%INTERP_HERMITE_LIM 1D limited Hermite interpolation
%
% [Vq,Dq] = INTERP_HERMITE_LIM(X,V,Xq,d_type) evaluates the limited Hermite
% interpolant Vq and its derivative Dq for a function V=F(X) at query
% points Xq using gradient estimate determined by the string d_type
% passed into calc_gradients (default 'hyman').
%
% [Vq,Dq] = INTERP_HERMITE_LIM(X,V,Xq) same as above with d_type='hyman'.
%
% X and V are vectors of the same length, and Vq and Dq will be returned
% with the same shape as Xq. For any points in Xq which are also in X, the
% corresponding Dq is the derivative from the right.
%
% See also: INTERP_HERMITE CALC_HERMITE EVAL_HERMITE

% Uses local subfunction apply_limiter.
%
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

D_lim = apply_limiter(X,V,D);

[Vq,Dq] = eval_hermite(X,V,D_lim,Xq);

end % function interp_hermite_lim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D_lim = apply_limiter(X,V,D)
%APPLY_LIMITER Restrict the derivatives for monotonicity
%
% Limits the derivatives for monotonicity, as described in Section 2.1 of
% Rasch and Williamson (1990), "On shape-preserving interpolation and
% semi-Lagrangian transport", SIAM J. Sci. Stat. Comput., Section 2.1. The
% resulting interpolating function is C0.

  grad = diff(V)./diff(X);

  D_lim = NaN(size(D));

  for ii = 1:length(X)-1
    D_l = D(ii,1);
    D_r = D(ii,2);

    % Test NCM0 for d(ii) from above
    if D(ii,1)*grad(ii) > 0
      % Impose SCM0
      D_l = sign(D_l)*min(abs(D_l),3*abs(grad(ii)));
    else
      D_l = 0;
    end % if
    % NCM0 for d(ii+1) from below
    if D(ii,2)*grad(ii) > 0
      % SCM0
      D_r = sign(D_r)*min(abs(D_r),3*abs(grad(ii)));
    else
      D_r = 0;
    end % if

    D_lim(ii,1) = D_l;
    D_lim(ii,2) = D_r;
  end % for ii
end % function apply_limiter
