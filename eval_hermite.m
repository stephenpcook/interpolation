function [Vq,Dq] = eval_hermite(X,V,D,Xq)
%EVAL_HERMITE Call a Hermite interpolant for each interval of X
%
% [Vq,Dq] = EVAL_HERMITE(X,V,D,Xq) returns the Hermite interpolant for data
% X, V and D, evaluated at query points Xq.
%
% Makes calls to subfunctions eval_hermite_local and eval_hermite_local_D
% for every interval of X.
%
% See also: INTERP_HERMITE INTERP_HERMITE_LIM

  Vq = NaN(size(Xq));
  Dq = NaN(size(Xq));

  N = length(X);
  for ii = 1:(N-1)
    x_l = X(ii);
    x_r = X(ii+1);
    v_l = V(ii);
    v_r = V(ii+1);

    % Not necessarially C1.
    d_l = D(ii,1);
    d_r = D(ii,2);

    % Just evaluate at the points in [x_l, x_r].
    % If Xq(jj) = X(ii+1) then it will be evaluated twice. Vq(jj) will remain
    % unchanged but Dq(jj) will take the value from the higher interval.
    local_mask = (x_l <= Xq) & (Xq <= x_r);
    Xq_loc = Xq(local_mask);

    Vq_loc = eval_hermite_local(x_l,x_r,v_l,v_r,d_l,d_r,Xq_loc);
    Dq_loc = eval_hermite_local_D(x_l,x_r,v_l,v_r,d_l,d_r,Xq_loc);

    Vq(local_mask) = Vq_loc;
    Dq(local_mask) = Dq_loc;
  end % for ii
end % function eval_hermite

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Vq = eval_hermite_local(x_l,x_r,v_l,v_r,d_l,d_r,Xq)
%EVAL_HERMITE_LOCAL 1D Hermite interpolant for one interval

% Hermite basis functions

  H1 = (x_r - Xq).^2/(x_r - x_l).^2 ...
       + (2*(Xq - x_l).*(x_r - Xq).^2)/((x_r - x_l)^3);
  H2 = (Xq - x_l).^2/(x_r - x_l).^2 ...
       + (2*(x_r - Xq).*(Xq - x_l).^2)/((x_r - x_l)^3);
  H3 = (Xq - x_l).*(x_r - Xq).^2/((x_r - x_l)^2);
  H4 = -(Xq - x_l).^2.*(x_r - Xq)/((x_r - x_l)^2);

  Vq = v_l*H1 + v_r*H2 + d_l*H3 + d_r*H4;

end % function eval_hermite_local

function Dq = eval_hermite_local_D(x_l,x_r,v_l,v_r,d_l,d_r,Xq)
%EVAL_HERMITE_LOCAL_D Derivative of the Hermite interpolant, one interval

  % Derivative of the basis functions from eval_hermite_local.
  H1_D = - (2*x_r - 2*Xq)/(x_l - x_r)^2 ...
         - (2*(x_r - Xq).^2)/(x_l - x_r)^3 ...
         - ((2*x_l - 2*Xq).*(2*x_r - 2*Xq))/(x_l - x_r)^3;
  H2_D = (2*(x_l - Xq).^2)/(x_l - x_r)^3 ...
         - (2*x_l - 2*Xq)/(x_l - x_r)^2 ...
         + (2*(x_r - Xq).*(2*x_l - 2*Xq))/(x_l - x_r)^3;
  H3_D = (x_r - Xq).^2/(x_l - x_r)^2 ...
         + ((x_l - Xq).*(2*x_r - 2*Xq))/(x_l - x_r)^2;
  H4_D = (x_l - Xq).^2/(x_l - x_r)^2 ...
         + ((x_r - Xq).*(2*x_l - 2*Xq))/(x_l - x_r)^2;

  Dq = v_l*H1_D + v_r*H2_D + d_l*H3_D + d_r*H4_D;

end % function eval_hermite_local_D
