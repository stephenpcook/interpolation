function mout = interp3lim(x, M, y)

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
