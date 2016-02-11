function pp = interp_ENO(X,F)
%INTERP_ENO Calculates the cubic ENO interpolant
%
% pp = INTERP_ENO(X,F) returns a piecewise polynomial structure of the
% cubic ENO interpolant with f(X(i)) = F(i). The coefficients for each
% interval are calculated with the function interp_ENO_coefs.
%
% See also: INTERP_ENO_COEFS INTERP1CUBICL INTERP3LIM PPVAL_LIM

% Stephen P. Cook 11-02-2016

pp_coefs = zeros(length(X)-1,4);

for ii = 1:length(X)-1
pp_coefs(ii,:) = interp_ENO_coefs(F,X,ii);
end

pp = mkpp(X, pp_coefs);

end % function interp_ENO
