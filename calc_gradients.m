function D = calc_gradients(X,V,d_type)
%CALC_GRADIENTS Wrapper function for calculating gradients
%
% D = CALC_GRADIENTS(X,V,d_type) numerically calculates the gradients the
% the data with the given method d_type. If X and V are of length N, then
% this will return the (N-1)-by-2 array with D(i,1) and D(i,2) are the left
% and right gradients of the interval [X(i),X(i+1)].
%
% d_type can be of the following:
%
%   'hyman'     - Quartic based on 4 surrounding points. From Hyman 1983
%   'akima'     - From Akima 1970
%   'quadratic' - Quadratic estimate based on 2 surrounding points
%   'zeros'     - All equal to zero

switch lower(d_type)
  case('hyman')
    D = hyman_gradients(X,V);
  case('akima')
    D = akima_gradients(X,V);
  case('quadratic')
    D = quadratic_gradients(X,V);
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

end % function calc_gradients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = quadratic_gradients(X,V)
%QUADRATIC_GRADIENTS Quadratic derivative esimates
%
% D = QUADRATIC_GRADIENTS(X,V) estimates the gradient of v(x) at x_i with a
% three point stencil.
%
% See also: AKIMA_GRADIENTS HYMAN_GRADIENTS

N = length(X);

Dx = diff(X);
S = diff(V)./Dx;

d(1) = ((2*Dx(1) + Dx(2))*S(1) - Dx(1)*S(2)) ...
         /(Dx(1)+Dx(2));
d(N) = ((2*Dx(N-1) + Dx(N-2))*S(N-1) - Dx(N-1)*S(N-2)) ...
             /(Dx(N-2)+Dx(N-1));
for ii = 2:N-1
  d(ii) = ( -Dx(ii)^2*V(ii-1) + (Dx(ii)^2-Dx(ii-1)^2)*V(ii) ...
            + Dx(ii-1)^2*V(ii+1)) ...
          ./(Dx(ii)*Dx(ii-1)*(Dx(ii) + Dx(ii-1)));
end
D = NaN(N-1,2);
D(:,1) = d(1:end-1);
D(:,2) = d(2:end);

end % function quadratic_gradients

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
% As described in Hyman (1983). Seems wrong at the moment

  N = length(X);

  grad = diff(V)./diff(X);

  % Special treatment near boundaries. "Second-order uncentred parabolic
  % method" from Hyman83 3C.
  d(1) = ((2*(X(2)-X(1)) + (X(3)-X(2)))*grad(1) - (X(2)-X(1))*grad(2)) ...
         /(X(3)-X(1));
  d(2) = ((2*(X(3)-X(2)) + (X(4)-X(3)))*grad(2) - (X(3)-X(2))*grad(3)) ...
         /(X(4)-X(2));
  for ii=3:N-2
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
              - (X(N-1)-X(N-2))*grad(N-3)) ...
             /(X(N-1) - X(N-3));
    d(N) = ((2*(X(N)-X(N-1)) + X(N-1)-X(N-2))*grad(N-1) ...
              - (X(N)-X(N-1))*grad(N-2)) ...
             /(X(N) - X(N-2));
  end % for ii

  D = NaN(N-1,2);
  D(:,1) = d(1:end-1);
  D(:,2) = d(2:end);

end % function akima_gradients
