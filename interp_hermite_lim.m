function [Vq,Dq] = interp_hermite_lim(X,V,Xq,d_type)
%INTERP_HERMITE_LIM 1D limited Hermite interpolation
%
% [Vq,Dq] = INTERP_HERMITE_LIM(X,V,Xq,d_type) evaluates the limited Hermite
% interpolant Vq and its derivative Dq for a function V=F(X) at query
% points Xq using gradient estimate determined by the string d_type
% ('hyman' (default), 'akima' or 'zeros').
%
% [Vq,Dq] = INTERP_HERMITE_LIM(X,V,Xq) same as above with d_type='hyman'.
%
% X and V are vectors of the same length, and Vq and Dq will be returned
% with the same shape as Xq. For any points in Xq which are also in X, the
% corresponding Dq is the derivative from the right.

% Uses local subfunctions hyman_gradients, akima_gradients, apply_limiter,
% eval_hermite, eval_hermite_local and eval_hermite_local_D.

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

switch lower(d_type)
  case('hyman')
    D = hyman_gradients(X,V);
  case('akima')
    D = akima_gradients(X,V);
  case('test')
    D = zeros(length(X),2);
    D(3,2) = 1;
    D(4,1) = 1;
  case('zeros')
    D = zeros(length(X),2);
  otherwise
    warning(['Unknown derivative type.', ...
             'Using default (hyman).']);
    D = hyman_gradients(X,V);
end % switch d_type

D_lim = apply_limiter(X,V,D);

[Vq,Dq] = eval_hermite(X,V,D_lim,Xq);

end % function interp_hermite_lim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = hyman_gradients(X,V)
%HYMAN_GRADIENTS Hyman finite difference derivative estimates
%
% As described in Hyman (1983) as "Fourth order finite difference".

  N = length(X);
  d = NaN(N,1);

  % Third order uncentred finite difference for boundaries, Hyman83 3C
  d(1) = (-22*V(1) + 36*V(2) - 18*V(3) + 4*V(4)) ...
           ./(-22*X(1) + 36*X(2) - 18*X(3) + 4*X(4));
  d(2) = (-2*V(1) - 3*V(2) + 6*V(3) - V(4)) ...
           ./(-2*X(1) - 3*X(2) + 6*X(3) - X(4));
  for ii = 3:N-2
    % Fourth order finite difference approximation, Hyman83 Table 1.
    % Is fourth order for "smooth-enough" meshes.
    d(ii) = (V(ii-2) - 8*V(ii-1) + 8*V(ii+1) - V(ii+2)) ...
              ./(X(ii-2) - 8*X(ii-1) + 8*X(ii+1) - X(ii+2));
  end % for ii
  % As for d(1) and d(2)
  d(N-1) = (V(N-3) - 6*V(N-2) + 3*V(N-1) + 2*V(N)) ...
           ./(X(N-3) - 6*X(N-2) + 3*X(N-1) + 2*X(N));
  d(N) = (-4*V(N-3) + 18*V(N-2) - 36*V(N-1) + 22*V(N)) ...
         ./(-4*X(N-3) + 18*X(N-2) - 36*X(N-1) + 22*X(N));

  D = NaN(N-1,2);
  D(:,1) = d(1:end-1);
  D(:,2) = d(2:end);
end % function hyman_gradients

function D = akima_gradients(X,V)
%AKIMA_GRADIENTS Akima derivative estimates
%
% As described in Hyman (1983).

  N = length(X);

  grad = diff(V)./diff(X);

  % Special treatment near boundaries. "Second-order uncentred parabolic
  % method" from Hyman83 3C.
  d(1) = ((2*(X(2)-X(1)) + (X(3)-X(2)))*grad(1) - (X(2)-X(1))*grad(2)) ...
         /(X(3)-X(1));
  d(2) = ((2*(X(3)-X(2)) + (X(4)-X(3)))*grad(2) - (X(3)-X(2))*grad(3)) ...
         /(X(4)-X(2));
  for ii=3:N-3
    % Akima estimate: Hyman83 Table 1 and RaschWilliamson90 Table 1.
    alpha_grad = abs(grad(ii+1) - grad(ii));
    beta_grad = abs(grad(ii-1) - grad(ii-2));
    if (alpha_grad + beta_grad == 0)
      d(ii) = grad(ii-1)+grad(ii);
    else
      d(ii) = (alpha_grad*grad(ii-1) + beta_grad*grad(ii)) ...
              /(alpha_grad + beta_grad);
    end % if
    d(N-1) = ((2*(X(N-1)-X(N-2)) + X(N-2)-X(N-3))*grad(N-2) ...
              - (X(N-2)-X(N-1))*grad(N-3)) ...
             /(X(N-1) - X(N-3));
    d(N) = ((2*(X(N)-X(N-1)) + X(N-1)-X(N-2))*grad(N-1) ...
              - (X(N-1)-X(N))*grad(N-2)) ...
             /(X(N) - X(N-2));
  end % for ii

  D = NaN(N-1,2);
  D(:,1) = d(1:end-1);
  D(:,2) = d(2:end);

end % function akima_gradients

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

function [Vq,Dq] = eval_hermite(X,V,D,Xq)
%EVAL_HERMITE Call a Hermite interpolant for each interval of X
%
% [Vq,Dq] = EVAL_HERMITE(X,V,D,Xq) returns the Hermite interpolant for data
% X, V and D, evaluated at query points Xq.
%
% Makes calls to eval_hermite_local and eval_hermite_local_D for every
% interval of X.

  Vq = NaN(size(Xq));
  Dq = NaN(size(Xq));

  N = length(X);
  for ii = 1:(N-1)
    x_l = X(ii);
    x_r = X(ii+1);
    v_l = V(ii);
    v_r = V(ii+1);

    % Not necessarially C1.
    d_l = D(ii,1);
    d_r = D(ii,2);

    % Just evaluate at the points in [x_l, x_r].
    % If Xq(jj) = X(ii+1) then it will be evaluated twice. Vq(jj) will remain
    % unchanged but Dq(jj) will take the value from the higher interval.
    local_mask = (x_l <= Xq) & (Xq <= x_r);
    Xq_loc = Xq(local_mask);

    Vq_loc = eval_hermite_local(x_l,x_r,v_l,v_r,d_l,d_r,Xq_loc);
    Dq_loc = eval_hermite_local_D(x_l,x_r,v_l,v_r,d_l,d_r,Xq_loc);

    Vq(local_mask) = Vq_loc;
    Dq(local_mask) = Dq_loc;
  end % for ii
end % function eval_hermite

function Vq = eval_hermite_local(x_l,x_r,v_l,v_r,d_l,d_r,Xq)
%EVAL_HERMITE_LOCAL 1D Hermite interpolant for one interval

% Hermite basis functions

  H1 = (x_r - Xq).^2/(x_r - x_l).^2 ...
       + (2*(Xq - x_l).*(x_r - Xq).^2)/((x_r - x_l)^3);
  H2 = (Xq - x_l).^2/(x_r - x_l).^2 ...
       + (2*(x_r - Xq).*(Xq - x_l).^2)/((x_r - x_l)^3);
  H3 = (Xq - x_l).*(x_r - Xq).^2/((x_r - x_l)^2);
  H4 = -(Xq - x_l).^2.*(x_r - Xq)/((x_r - x_l)^2);

  Vq = v_l*H1 + v_r*H2 + d_l*H3 + d_r*H4;

end % function eval_hermite_local

function Dq = eval_hermite_local_D(x_l,x_r,v_l,v_r,d_l,d_r,Xq)
%EVAL_HERMITE_LOCAL_D Derivative of the Hermite interpolant, one interval

  % Derivative of the basis functions from eval_hermite_local.
  H1_D = - (2*x_r - 2*Xq)/(x_l - x_r)^2 ...
         - (2*(x_r - Xq).^2)/(x_l - x_r)^3 ...
         - ((2*x_l - 2*Xq).*(2*x_r - 2*Xq))/(x_l - x_r)^3;
  H2_D = (2*(x_l - Xq).^2)/(x_l - x_r)^3 ...
         - (2*x_l - 2*Xq)/(x_l - x_r)^2 ...
         + (2*(x_r - Xq).*(2*x_l - 2*Xq))/(x_l - x_r)^3;
  H3_D = (x_r - Xq).^2/(x_l - x_r)^2 ...
         + ((x_l - Xq).*(2*x_r - 2*Xq))/(x_l - x_r)^2;
  H4_D = (x_l - Xq).^2/(x_l - x_r)^2 ...
         + ((x_r - Xq).*(2*x_l - 2*Xq))/(x_l - x_r)^2;

  Dq = v_l*H1_D + v_r*H2_D + d_l*H3_D + d_r*H4_D;

end % function eval_hermite_local_D
