function ts = calculateSettlementTime(e, t, tol) 
    %%%Eliminate NAN values
    e(isnan(e)) = e(end-1,:);
    %%%Set bounds
    bounds = [-tol tol];
    %%%find index out the frange
    [x,~] = find(~(e<bounds(2) & e>bounds(1)));
    if(x(end)==length(t))
        ts=NaN;
    else
        %%%find the last index out the frange
        ts = t(max(x));
    end
end