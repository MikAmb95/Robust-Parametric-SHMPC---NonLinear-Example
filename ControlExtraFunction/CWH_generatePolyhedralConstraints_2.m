function [xLim, uLim] =  CWH_generatePolyhedralConstraints_2()
nu = 3;
nx = 6;


% input constraints
uLim.b = ones(nu*2,1)*0.5; %m/s^2
uLim.A = [-eye(nu); eye(nu)];

%state constraints

% min max values for state box constraints
xBox =  ([[-ones(3,1)*10000 ones(3,1)*10000]; ...
    [-ones(3,1)*10000 ones(3,1)*10000]]);    % adding some constraints to have a compact set


%linear constraints derived from saturation values / box
xLim.A =kron([-1;1],eye(nx));
xLim.b= [-xBox(:,1); xBox(:,2)];

end