function [rho_pp, rho_int_pp] = slice_pp2(z_A, zh_A, rhobar)

% global N
N = 60;
%=======
% INPUT
%=======
% z_A
% zh_A
% rhobar
xi = (0:N)./N;
z_ceiling = 10000;
alpha_Q = 0.3;
z_A = (alpha_Q*xi.^2 + (1-alpha_Q)*xi) * z_ceiling;
zh_A = (z_A(1:end-1) + z_A(2:end))/2;
% rhobar with a perturbation
rhobar = ones(1,N)-5e-5*zh_A.*cos(2*pi*zh_A/z_ceiling);
zmid_i = floor(N/2) + (0:3);
rhobar(zmid_i) = 1.3*rhobar(zmid_i);


%========
% OUTPUT
%========
% rho_pp
% rho_int_pp

% Initialisation
%----------------
% Number of intervals of the pp
Dz_A = diff(z_A);
intervals = N-4;


% Okay, for each interval that we work in, we want to create the pp_coefs
% of the mass functions with Lagrangian polynomials.

masses = Dz_A.*rhobar;
% coefs of a lagrangian polynomial?

% We need to evaluate rho_L and rho_R (derivative of cumulative masses) at
% all x_i's and at all midpoints, and rho_L_prime and rho_R_prime (second
% derivative of cumulative masses) at all the midpoints.

%%%%%%%%%%%%%
% MAIN LOOP %
%%%%%%%%%%%%%

for ii = 2:(n-2)
  % All assuming indexing starts at zero for both X and masses
  X_local = X(ii-3):X(ii+1);
  Y_local = [0; cumsum(masses(ii-2:ii+1))];
  X_halfs = [(X(ii-2)+X(ii-1))/2, (X(ii-1)+X(ii))/2];

  rho_L(ii) = diff_lagrange(X_local, Y_local, X(ii-1))
  [rho_L_prime_m(ii), rho_L_prime_p(ii)] = ...
    ddiff_quartic_lagrange(X_local,Y_local,X_halfs);
end % for ii
% This should be the code for calculating the defining points
rho_i_minus = rho_L(i-1);
rho_i = rho_L(i);
rho_i_half = 3*masses(i)/(2*Dxi(i)) - 0.25*(rho_i_minus + rho_i);
rho_i_plus = calc_rho_i_plus_one(rho_i_minus, rho_i_half, rho_i,...
  rho_L_prime_p(i), rho_L_prime_m(i+1), Dxi, Dxi_plus);
% These above should be calculated once then output (quite expensive), and
% then they can be used many times (relatively cheap).
%
%rho_D = my_lagrange(%[X(i-1), 0.5*(X(i-1)+X(i)), X(i), X(i+1)],...
%                    %[rho_i_minus, rho_i_half, rho_i, rho_i_plus]...
%                    %,X_out)

function p_three = calc_rho_i_plus_one(rho_i_minus, rho_i_half, rho_i,...
                                d_rho_L_half, d_rho_R_half, Dxi, Dxi_plus)
  % Calculates p_three as in thesis/slice.tex equation eq:slice:pThree
  pi_prime_i_m = (0.5*Dxi + Dxi_plus)/(Dxi*(Dxi + Dxi_plus));
  pi_prime_i_half = 1/(0.5*Dxi + Dxi_plus);
  pi_prime_i = (0.5*Dxi + Dxi_plus)/(Dxi*Dxi * plus);
  pi_prime_i_plus = 1/4*(Dxi^2)/...
    ((Dxi + Dxi_plus)*(0.5*Dxi+Dxi_plus)*Dxi_plus);

  p_three = 1/pi_prime_i_plus * (0.5*(d_rho_L_half + d_rho_R_half) - ...
    rho_i_minus*pi_prime_i_m - rho_i_half*pi_prime_i_half - ...
  rho_i*pi_prime_i);
end % function p_three
end % function slice_pp2
