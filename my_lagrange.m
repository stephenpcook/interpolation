function [Y_out] = my_lagrange(X_in,Y_in,X_out)
% Evaluates a Lagrangian polynomial for given values
%
%  X_in  - Points at which the polynomial is defined
%  Y_in  - Value of the function to be interpolated
%  X_out - Points at which the lagrange polynomial is to be evaluated
%
% [Y_out] = my_lagrange(X_in,Y_in,X_out)

%   $$f(x) = \sum_1^5 f_i * pi_ii(x)$$
% where
%   $$pi_ii = \prod_{j = 1 \\ j \ne i}^n (x-X_j)/(X_i - X_j)$$

n = length(X_in);

Y_out = zeros(size(X_out));
for ii = 1:n
  pi_ii = 1;
  for jj = 1:n
    if jj~=ii
      pi_ii = pi_ii.*(X_out - X_in(jj))./(X_in(ii) - X_in(jj));
    end % if jj ~= ii
  end % for jj
  Y_out = Y_out + Y_in(ii)*pi_ii;
end % for ii

end % function diff_lagrange
