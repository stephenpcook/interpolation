function YY = ppval_lim(pp, xx)

% Now evaluate the y's where they're supposed to be!
xx_N = length(xx);

YY= ppval(pp,xx);

[BREAKS,COEFS,L,K,D] = unmkpp(pp);
% The value at the right-most break. Only one not included in 
% the pp format.
y_end = polyval(COEFS(end,:),BREAKS(end)-BREAKS(end-1));


for jj = 1:xx_N
    % Find the correct interval. Could also be done with histc?
    interval = 1; 
    while xx(jj)>BREAKS(interval + 1) 
        interval = interval + 1;
    end % while
    y_l = COEFS(interval,K);
    if interval>(L-1) 
        % Okay, this is a slight hack. We done have the value at the
        % right-most break.
        y_r = y_end;
    else
        y_r = COEFS(interval+1,K);
    end
    y_max = max(y_l,y_r);
    y_min = min(y_l,y_r);
    if YY(jj) > y_max
        YY(jj) = y_max;
        % Warning messages if appropriate verbatim param
        %warning(['Overshoot limited in interval (',num2str(x(interval)),...
        %    ', ',num2str(x(interval+1)),')'])
    elseif YY(jj) < y_min
        YY(jj) = y_min;
        %warning(['Undershoot limited in interval (',num2str(x(interval)),...
        %   ', ',num2str(x(interval+1)),')'])
        % Warning messages if appropriate verbatim param
    end % if
end % for jj=1:y_N
end % function