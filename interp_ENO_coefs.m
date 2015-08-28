function pp_coefs = interp_ENO_coefs(F,X,i_start)
% Calculates the cubic lagrange ENO interpolating function coefficients
%
% For the interval [X(i_start), X(i_start + 1)], this finds the escentially
% non-oscillatory cubic (ENO) lagrange interpolating function for F defined
% on either X([i-2,i-1,i,i+1]), X([i-1,i,i+1,i+2]), or X([i,i+1,i+2,i+3])
% which minimises oscillations [HartenEtAl86, Smith00].
%   The function is found via the divided difference formulation, then the
% coefficients for a MATLAB piecewise polynomial are found.
%   Naming convention is
%     F1N := f[X_1,X_2,...,X_N]
%          = (f[X_2,...,X_N] - f[X_1,...,X_N-1])/(X_N-X_1)
% and
%      FN := f[X_N] = f(X_N).
%
% IN
%   F      - Function values to be interpolated
%   X      - Points corresponding to F
%  i_start - Left side of interval over which to interpolate.
%
% OUT
%   pp_coefs - [c1,c2,c3,c4], coefficients of the interpolating cubic.
%
% Subfunctions
%   second_diff
%   third_diff
%   poly_coefs
%
% Stephen Cook 14/05/2015
% 
% References
%   [HartenEtAl86] A. Harten, S. Osher, B. Engquist, and S.R. Chakravarthy. 
%           Some results on uniformly high order accurate  essentially
%           non-oscillatory schemes. Applied. Num. Math., 2:347-377, 1986.
%   [Smith00] C.J. Smith, The semi-Lagrangian method in atmospheric 
%           modelling, Ph.D. thesis, University of Reading, 2000.

%%%%%%%%%%%%%%%%%%%%%
% Main body of code %
%%%%%%%%%%%%%%%%%%%%%

% Setup
% i_bot and i_top mark the highest and lowest indices of X in X_
% X_ and F_ are the values of X and F currently used for divided
% difference.
i_bot = i_start;
i_top = i_start + 1;
X_ = [X(i_bot),X(i_top)];
F_ = [F(i_bot),F(i_top)];

% Zero-th divided differences.
F1 = F_(1);
F2 = F_(2);

% First divided difference.
F12 = (F2 - F1)/(X_(2) - X_(1));

% Second divided difference, with ENO test.
% Testing the point to either side to see which gives the least oscilation.
% Note that X_ will not necessarially be ordered, which isn't an issue for
% divided difference.
%
% Only one direction to go if we're at the top or bottom of X.
if i_bot==1
    i_top = i_top + 1;
    X_ = [X_,X(i_top)];
    F_ = [F_,F(i_top)];
    [F3,F23,F123] = second_diff(F2,F12,F_,X_);
elseif i_top==length(X)
    i_bot = i_bot - 1;
    X_ = [X_,X(i_bot)];
    F_ = [F_,F(i_bot)];
    [F3,F23,F123] = second_diff(F2,F12,F_,X_);
else 
    % Test point above and point below
    Xp = X(i_top+1);
    Fp = F(i_top+1);
    [F3p,F23p,F123p] = second_diff(F2,F12,[F_,Fp],[X_,Xp]);
    Xm = X(i_bot-1);
    Fm = F(i_bot-1);
    [F3m,F23m,F123m] = second_diff(F2,F12,[F_,Fm],[X_,Xm]);
    if abs(F123m) <= abs(F123p)
        % Looks like we default to backwards difference if F123m==F123p.
        % Not an issue, as the fourth divided difference would give
        % F1234p==0.
        i_bot = i_bot - 1;
        X_ = [X_,Xm];
        F_ = [F_,Fm];
        F3 = F3m;
        F23 = F23m;
        F123 = F123m;
    else
        i_top = i_top + 1;
        X_ = [X_,Xp];
        F_ = [F_,Xp];
        F3 = F3p;
        F23 = F23p;
        F123 = F123p;
    end
end


% Third divided difference with ENO test.
% Since it's the last divided difference for cubics, don't actually need
% all the values; commented out.
if i_bot==1
    i_top = i_top + 1;
    X_ = [X_,X(i_top)];
    F_ = [F_,F(i_top)];
    %[F4,F34,F234,F1234] = third_diff(F3,F23,F123,F_,X_);
    [~,~,~,F1234] = third_diff(F3,F23,F123,F_,X_);
elseif i_top==length(X)
    i_bot = i_bot - 1;
    X_ = [X_,X(i_bot)];
    F_ = [F_,F(i_bot)];
    %[F4,F34,F234,F1234] = third_diff(F3,F23,F123,F_,X_);
    [~,~,~,F1234] = third_diff(F3,F23,F123,F_,X_);
else
    Xp = X(i_top+1);
    Fp = F(i_top+1);
    %[F4p,F34p,F234p,F1234p] = third_diff(F3,F23,F123,[F_,Fp],[X_,Xp]);
    [~,~,~,F1234p] = third_diff(F3,F23,F123,[F_,Fp],[X_,Xp]);
    Xm = X(i_bot-1);
    Fm = F(i_bot-1);
    %[F4m,F34m,F234m,F1234m] = third_diff(F3,F23,F123,[F_,Fm],[X_,Xm]);
    [~,~,~,F1234m] = third_diff(F3,F23,F123,[F_,Fm],[X_,Xm]);
    if abs(F1234m) <= abs(F1234p)
        %i_bot = i_bot - 1;
        X_ = [X_,Xm];
        %F_ = [F_,Fm];
        %F4 = F4m;
        %F34 = F34m;
        %F234 = F234m;
        F1234 = F1234m;
    else
        %i_top = i_top + 1;
        X_ = [X_,Xp];
        %F_ = [F_,Xp];
        %F4 = F4p;
        %F34 = F34p;
        %F234 = F234p;
        F1234 = F1234p;
    end % if abs(F1234m) <= abs(F1234p)
end

% Output

% Coefficients for mkpp relative to X(i_start),
pp_coefs = poly_coefs(F1,F12,F123,F1234,X_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of main body of code %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% Subfunctions %
%%%%%%%%%%%%%%%%

function [F3,F23,F123] = second_diff(F2,F12,F_,X_)
  % Calculate the second divided difference
  %   F123 := f[x1,x2,x3]
  %         = (f[x2,x3] - f[x1,x2])/(x3-x1) .
    F3 = F_(3);
  F23 = (F3 - F2)/(X_(3) - X_(2));
  F123 = (F23 - F12)/(X_(3) - X_(1));
end % subfunction second_diff

function [F4,F34,F234,F1234] = third_diff(F3,F23,F123,F_,X_)
  % Calculate the third divided difference
  %   F1234 := f[x1,x2,x3,x4]
  %          = (f[x2,x3,x4] - f[x1,x2,x3])/(x4-x1) .
  F4 = F_(4);
  F34 = (F4 - F3)/(X_(4) - X_(3));
  F234 = (F34 - F23)/(X_(4) - X_(2));
  F1234 = (F234 - F123)/(X_(4) - X_(1));
end % subfunction third_diff

function [out] = poly_coefs(F1,F12,F123,F1234,X_)
% Gives the coefficient of the polynomial
%   p3(x) = c1*(x-X_1)^3 + c2*(x-X_1)^2 + c3*(x-X_1) + c4 .
% These coefficients were arrived at by expanding the divided difference
% equation,
%   p3(x) = F1 + (x-X1)*F12 + (x-X1)*(x-X2)*F123 + (x-X1)*(x-X2)*(x-X3)*F1234
% which, in a local coordinate system (relative to X1) is
%   p3(x) = F1 + x*F12 + x*(x-d12)*F123 + x*(x-d12)*(x-d13)*F1234
% with 
%   d12 = X2 - X1,
%   d13 = X3 - X1.
d12 = X_(2) - X_(1);
d13 = X_(3) - X_(1);
    
c4 = F1;
c3 = F12 - d12*(F123 - d13*F1234);
c2 = F123 - (d12 + d13)*F1234;
c1 = F1234;
out = [c1,c2,c3,c4];
end % subfunction poly_coefs
end % function
