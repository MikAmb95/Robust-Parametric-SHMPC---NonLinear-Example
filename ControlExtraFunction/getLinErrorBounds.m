function [LE_bounds,LE_bounds_interval,gamma] = getLinErrorBounds(nx,Z,zstar)

    intZ = interval(Z);

    gamma = max(abs(zstar-intZ.inf),abs(intZ.sup-zstar));
    f = my_H_handle(nx);


    L = zeros(nx,1);
    for i=1:nx
        hfunc = f{i};
        %evaluates hessian over interval set and computes absolute value
        %the outer interval conversion is there just for cases where the
        %hfunc function evaluation returns a matrix of doubles even though the input is an interval, for
        %example, if all elements of hfunc are zero.
        hessInt = abs(interval(hfunc(intZ)));
        hessMax = hessInt.sup; %maximum absolute value of the hessian over DeltaZ
        L(i)  = 0.5*gamma'*hessMax*gamma;
    end

    LE_bounds = zonotope(zeros(nx,1),diag(L));
    LE_bounds_interval= interval(LE_bounds);
end