function mout = interp3lim(x, M, y)
%INTERP3LIM cubic lagrange imterpolation with a flux limiter
%
% INTERP3LIM(x, M, y) uses interp1cubicL then uses a flux limiter such that
% the interpolant is piecewise monotonic.
%
% Example: As in the interp1cubicL example, we can compare the two
% functions (interp1cubicL blue dashed, interp3lim red solid) with
%
%     X = (0:6);
%     M = sin(X);
%     y = 0:0.1:6;
%     plot(y,interp1cubicL(X,M,y),'b--', y,interp3lim(X,M,y),'r-', X,M,'b+')
%
% See also: INTERP1CUBICL PPVAL_LIM PCHIP

% How to treat extrapolation?
% Now evaluate the y's where they're supposed to be!
mout= interp1cubicL(x,M,y);
y_N = length(y);
for jj = 1:y_N
    interval = 1;
    if y(jj)>x(end)
        interval = length(M)-1;
        %continue
    else
    % Find the right interval
    while y(jj)>x(interval + 1)
        interval = interval + 1;
    end % while
    end %if
    Mmax = max(M(interval),M(interval+1));
    Mmin = min(M(interval),M(interval+1));
    if mout(jj) > Mmax
        mout(jj) = Mmax;
        % Warning messages if appropriate verbatim param
        %warning(['Overshoot limited in interval (',num2str(x(interval)),...
        %    ', ',num2str(x(interval+1)),')'])
    elseif mout(jj) < Mmin
        mout(jj) = Mmin;
        %warning(['Undershoot limited in interval (',num2str(x(interval)),...
        %   ', ',num2str(x(interval+1)),')'])
        % Warning messages if appropriate verbatim param
    end % if
end % for jj=1:y_N
end % function
