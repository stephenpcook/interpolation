% function [rho_pp, rho_int_pp] = slice_pp2(z_A, zh_A, rhobar)

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
