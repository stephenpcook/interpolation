function [Y_prime_out] = diff_lagrange(X_in,Y_in,X_out)
% Calculates derivatives of a Lagrangian polynomial for given values
%
%  X_in  - Points at which the polynomial is defined
%  Y_in  - Value of the function to be interpolated
%  X_out - Points at which the derivative is to be evaluated
%
% [Y_prime_out] = diff_lagrange(X_in,Y_in,X_out)


%f(x) = sum_1^n f_i * pi_ii(x);
%f_prime = sum_1^n f_i * pi_i_prime(x);
%
% pi_ii = prod_{j = 1 \\ j \ne i}^n (x-X_j)/(X_i - X_j)
% pi_ii_prime = \sum_{j=1 \\ j \ne i}^n
%                \prod_{k=1 \\ k \ne j \\ k \ne i}^n (x-X_k)/(X_i - X_k)

n = length(X_in);

Y_prime_out = zeros(size(X_out));
for ii = 1:n;
  pi_ii_prime_numerator = zeros(size(X_out));
  for kk = 1:n
    % Summation loop
    if kk ~= ii
      pi_ii_kk = 1;
      for jj = 1:n
        % Product loop
        if (jj ~= ii) & (jj ~= kk)
          pi_ii_kk = pi_ii_kk.*(X_out - X_in(jj));%./(X_in(ii) - X_in(jj));
        end % if jj ~= ii
      end % for jj
      pi_ii_prime_numerator = pi_ii_prime_numerator + pi_ii_kk;
    end % if kk ~= ii
  end % for kk
  pi_ii_prime_denominator = 1;
  for jj = 1:n
    if jj ~= ii
      pi_ii_prime_denominator = pi_ii_prime_denominator*(X_in(ii) - X_in(jj));
    end % if jj ~= ii
  end % for jj
  Y_prime_out = Y_prime_out + Y_in(ii)*pi_ii_prime_numerator/pi_ii_prime_denominator;
end % for ii

end % function diff_lagrange
