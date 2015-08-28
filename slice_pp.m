function [rho_pp, rho_int_pp] = slice_pp(z_A, zh_A, rhobar)
% Creates a mass preserving pp of rho and its integral.

% This should be part of a set of functions for performing the SLICE
% interpolation [Zerroukat Et Al. 2002], where-by a piecewise polynomial is
% constructed from the masses advected from the previous timestep, and a
% function for doing this advection.
%
% An extension would be to implement the monotone version of this from
% [Zerroukat Et Al. 2005]. There's room for MUCH improvement, but that
% might involve reformulating the definition of rhohat without the local
% transformation...
%
% TODO - remove the assumption that zh_A are midpoints of z_A?


global N
%=======
% INPUT 
%=======
% z_A - The points defining the eularian control volumes (ECV).
% zh_A - Currently not used; assumed to be the midpoints of z_A.
% rhobar - is the mass in each ECV. Sit at zh_A.


%========
% OUTPUT
%========
% rho_pp - The mass preserving piecewise polynomial (pp).
% rho_int_pp - The exact integral of rho_pp.

% Initialisation
%----------------
% Number of intervals of the pp
intervals = N-4;

hat_coefs = zeros(intervals+1, 4);
rho_coefs = zeros(intervals, 4);
int_coefs = zeros(intervals, 5);

Dz_A = diff(z_A);

% rho_hat
%---------
% First we want to build up the N-3 interval pp rhohat(z) with 
%   \int_{I_j} rhohat(z) = rhobar_j
% with j \in \{i-2, i-1, i, i+1\}.
%
%  Values used in the matrix.
C1 = z_A(1:end-1) + z_A(2:end);
C2 = z_A(1:end-1).^2 + z_A(2:end).^2;
C3 = z_A(1:end-1).*z_A(2:end);
C1 = C1(:); 
C2 = C2(:);
C3 = C3(:); 

% We're going to move through this pulling out a 4-by-4 matrix from hugeA
% for each interval. Calculated from analytical integrals required to meet
% the conditions above.
hugeA = [C1.*C2./4, (C2+C3)./3, C1./2, ones(N,1)];

% Local coefficients of the rho_hat pp's.
%-----------------------------------------
for ii = 1:N-3; % Create all the coefficients for the pp of rhohat.
  
  hat_coefs(ii,1:4) = (hugeA(ii:ii+3, 1:4)\rhobar(ii:ii+3)')';
  
  % What a horrible mess. So the coefficients I found were for a regular
  % polynomial, like
  %   if z_L < z < z_R
  %     pp1(z) = A1*z^3 + A2*z^2 + A3*z + A4
  %   end if
  % but MATLAB actually does it as a local polynomial,
  %   if z_L < z < z_R
  %     pp1(z) = a1*(z-z_L)^3 + a2*(z-z_L)^2 + a3*(z-z_L) + a4
  %   end if
  % Instead of going back to derive the a#'s corresponding to this
  % formulation, I'm just going to do a global to local transformation. 
  %  I'll call this an acceptable intellectual overhead.
  zi = z_A(ii+2);
  local_change = [...
    1, 0, 0, 0;...
    3*zi, 1, 0, 0;...
    3*zi^2, 2*zi, 1, 0;...
    zi^3, zi^2, zi, 1];
  
  % Row vector, so right hand multiplication
  hat_coefs(ii,1:4) = hat_coefs(ii,1:4)*local_change';
end

%hat_breaks = z_A(3:end-1);
%rhohat= mkpp(hat_breaks, hat_coefs);

% Making the pp of rho from the rho_hat pp.
%-------------------------------------------
% Now we're going to make the other piecewise polynomial, defined by
%   rho(z_{i-1/2}) = rhohat_i(z_{i-1/2}),
%   rho(z_{i+1/2}) = rhohat_{i+1}(z_{i+1/2}),
%   \int_{I_i} rho(z)dz = rhobar_i,
%   d rho/dz |_{z_i} = 1/2 (d rhohat_i/dz + d rhohat_{i+1}/dz)|_{z_i}.
for ii = 1:intervals
  dz = Dz_A(ii+2);
  A = [...
    0, 0, 0, 1;...
    dz^3, dz^2, dz, 1;...
    1/4*dz^4, 1/3*dz^3, 1/2*dz^2, dz;...
    3/4*dz^2, dz, 1, 0];
  b = [...
    hat_coefs(ii, 4);...
    hat_coefs(ii+1, 4);...
    rhobar(ii+2)*dz;...
    1/1*(3/8*dz^2*(hat_coefs(ii,1) + hat_coefs(ii+1,1)) + ...
      1/2*dz*(hat_coefs(ii,2) - hat_coefs(ii+1,2)) + ...
      1/2*(hat_coefs(ii,3) + hat_coefs(ii+1,3)))];
  rho_coefs(ii,1:4) = (A\b)';
end % for ii = 1:intervals
rho_breaks = z_A(3:end-2);

% Make the pp of the integral of rho.
%------------------------------------
int_coefs(:,1:4) = rho_coefs .* (ones(intervals,1) * [0.25, 1/3, 0.5, 1]);

% Fix the constants of integration to get the definite integral (continuous)
mass_sum = cumsum(rhobar(3:end-3).*diff(rho_breaks(2:end)));
int_coefs(:,5) = [0; mass_sum(:)];

% Construct the pps
%-------------------
rho_pp = mkpp(rho_breaks, rho_coefs);
rho_int_pp = mkpp(rho_breaks, int_coefs);

end % function slice_interp
