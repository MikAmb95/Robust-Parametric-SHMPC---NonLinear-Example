
function E = computeNStepE_lin(nx,Aerr,N,W,K,z_star)
k = 0;
vx0 = zeros(nx,1); 
options.zonotopeOrder = 6;
options.reductionTechnique = 'girard';
E{1} = zonotope(interval(vx0,vx0));
E_sum = zonotope(interval(vx0,vx0));
KE = mtimes(K,E{1});
term1 = cartProd(E{1},cartProd(KE,W));
Z{1} = z_star + term1;
[LE_bounds,~,~] = getLinErrorBounds(nx,Z{1},z_star);
for i = k+1:N-1
    E_sum =  mtimes(Aerr,E{i}) + W + LE_bounds;
    E_sum = reduce(E_sum,options.reductionTechnique,options.zonotopeOrder);
    E{i+1} = E_sum;
    KE = mtimes(K,E_sum);
    term1 = cartProd(E_sum,cartProd(KE,W));
    Z_sum = z_star + term1;
    Z_sum = reduce(Z_sum,options.reductionTechnique,options.zonotopeOrder);
    Z{i+1} = Z_sum;
    %[LE_bounds,~,~] = getLinErrorBounds(nx,Z{i+1},z_star);
end
E_sum = mtimes(Aerr,E{i}) + W;
E_sum = reduce(E_sum,options.reductionTechnique,options.zonotopeOrder);
E{N} = E_sum;
end

