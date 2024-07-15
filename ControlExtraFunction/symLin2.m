


function [A,B] = symLin2(x0, u0)

x = sym('x',[6,1],'real');
u = sym('u',[3,1],'real');

J1 = 120; J2 = 100; J3= 80;

phi = x(1); % roll
theta = x(2); % pitch
psi = x(3); % yaw

om1 = x(4); % first component of the angular velocity
om2 = x(5); % second component of the angular velocity
om3 = x(6); % third component of the angular velocity


% Kinematic equations
adot = 1/cos(theta)*[cos(theta), sin(phi)*sin(theta),cos(phi)*sin(theta);
0, cos(phi)*cos(theta), -sin(phi)*cos(theta);
0, sin(phi), cos(phi)]*[om1;om2;om3];


phidot = adot(1); % time rate of change of roll
thetadot = adot(2); % time rate of change of pitch
psidot = adot(3); % time rate of change of yaw

% External moments
M = u(:);
M1 = M(1);
M2 = M(2);
M3 = M(3);

% Dynamic equations
om1dot = ((J2-J3)*om2*om3)/J1+M1/J1;
om2dot = ((J3-J1)*om3*om1)/J2+M2/J2;
om3dot = ((J1-J2)*om1*om2)/J3+M3/J3;

dxdt = [phidot, thetadot, psidot, om1dot, om2dot, om3dot]';

A = jacobian(dxdt,x);
B = jacobian(dxdt,u);

A = subs(A,x,x0);
B = subs(B,u,u0);

end









