function pp = interp_ENO(X,F)

pp_coefs = zeros(length(X)-1,4);
%breaks = X;

for ii = 1:length(X)-1
%i_start = 2;
pp_coefs(ii,:) = interp_ENO_coefs(F,X,ii);
end

%%%%%%%%%%%%
% plotting %
%%%%%%%%%%%%
pp = mkpp(X, pp_coefs);

end % function interp_ENO