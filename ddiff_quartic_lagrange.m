function [Y_dd_out] = ddiff_quartic_lagrange(X_in,Y_in,X_out)
% Evaluates the second derivative of a quartic lagrange polynomial
%
%  X_in  - 5 Points at which the quartic is defined
%  Y_in  - 5 Value of the function to be interpolated
%  X_out - Points at which the second derivative is to be evaluated
%
% [Y_dd_out] = ddiff_quartic_lagrange(X_in,Y_in,X_out)

% Some initial checks
if length(X_in) ~= 5 & langth(Y_in) ~= 5
  error('Need 5 points to define a quartic');
end % if length(X_in)
if length(unique(X_in)) < 5
  error('X_in must be distinct');
end % if length(unique(X_in))

Y_dd_out = 0;
x = X_out;

X1 = X_in(1);
X2 = X_in(2);
X3 = X_in(3);
X4 = X_in(4);
X5 = X_in(5);

% Second derivatives of the lagrange interpolating functions
% (All been checked)
dd_pi_1 = 2*((x-X2).*(x-X3) + (x-X2).*(x-X4) + (x-X2).*(x-X5) + ...
             (x-X3).*(x-X4) + (x-X3).*(x-X5) + (x-X4).*(x-X5))./ ...
            ((X1-X2)*(X1-X3)*(X1-X4)*(X1-X5));
dd_pi_2 = 2*((x-X1).*(x-X3) + (x-X1).*(x-X4) + (x-X1).*(x-X5) + ...
             (x-X3).*(x-X4) + (x-X3).*(x-X5) + (x-X4).*(x-X5))./ ...
            ((X2-X1)*(X2-X3)*(X2-X4)*(X2-X5));
dd_pi_3 = 2*((x-X1).*(x-X2) + (x-X1).*(x-X4) + (x-X1).*(x-X5) + ...
             (x-X2).*(x-X4) + (x-X2).*(x-X5) + (x-X4).*(x-X5))./ ...
            ((X3-X1)*(X3-X2)*(X3-X4)*(X3-X5));
dd_pi_4 = 2*((x-X1).*(x-X2) + (x-X1).*(x-X3) + (x-X1).*(x-X5) + ...
             (x-X2).*(x-X3) + (x-X2).*(x-X5) + (x-X3).*(x-X5))./ ...
            ((X4-X1)*(X4-X2)*(X4-X3)*(X4-X5));
dd_pi_5 = 2*((x-X1).*(x-X2) + (x-X1).*(x-X3) + (x-X1).*(x-X4) + ...
             (x-X2).*(x-X3) + (x-X2).*(x-X4) + (x-X3).*(x-X4))./ ...
            ((X5-X1)*(X5-X2)*(X5-X3)*(X5-X4));

Y_dd_out = Y_in(1)*dd_pi_1 +  Y_in(2)*dd_pi_2 +  Y_in(3)*dd_pi_3 + ...
              Y_in(4)*dd_pi_4 +  Y_in(5)*dd_pi_5;


end % function ddiff_quartic_lagrange
