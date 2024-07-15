
function f = my_H_handle(nx)

%f{1} 
f{1} = @(x) [ -(sin(x(2))*(x(6)*cos(x(1)) + x(5)*sin(x(1))))/cos(x(2)),           (x(5)*cos(x(1)) - x(6)*sin(x(1)))/cos(x(2))^2, 0, 0, (cos(x(1))*sin(x(2)))/cos(x(2)), -(sin(x(1))*sin(x(2)))/cos(x(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                      x(6)*sin(x(1)) - x(5)*cos(x(1)),                                             0, 0, 0,                  -sin(x(1)),                   -cos(x(1)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
           -(x(6)*cos(x(1)) + x(5)*sin(x(1)))/cos(x(2)), (sin(x(2))*(x(5)*cos(x(1)) - x(6)*sin(x(1))))/cos(x(2))^2, 0, 0,           cos(x(1))/cos(x(2)),           -sin(x(1))/cos(x(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                                            0,                                             0, 0, 0,                         0,                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                                            0,                                             0, 0, 0,                         0,                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                                            0,                                             0, 0, 0,                         0,                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0;zeros(9,15)];
                                        


%%%%%%%%%

%f{2}

f{2} = @(x)[           (x(5)*cos(x(1)) - x(6)*sin(x(1)))/cos(x(2))^2,        (2*sin(x(2))*(x(6)*cos(x(1)) + x(5)*sin(x(1))))/cos(x(2))^3, 0, 0,           -sin(x(1))/(sin(x(2))^2 - 1),           cos(x(1))/cos(x(2))^2, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                                             0,                                                      0, 0, 0,                                  0,                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 (sin(x(2))*(x(5)*cos(x(1)) - x(6)*sin(x(1))))/cos(x(2))^2, -((cos(x(2))^2 - 2)*(x(6)*cos(x(1)) + x(5)*sin(x(1))))/cos(x(2))^3, 0, 0, -(sin(x(1))*sin(x(2)))/(sin(x(2))^2 - 1), (cos(x(1))*sin(x(2)))/cos(x(2))^2, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                                             0,                                                      0, 0, 0,                                  0,                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                                             0,                                                      0, 0, 0,                                  0,                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                                             0,                                                      0, 0, 0,                                  0,                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0;zeros(9,15)];
 




%%%%%%%%%%%%%%%

%f{3} 

f{3} = @(x) zeros(15,15);

%%%%%%%%%%%%%%%

%f{4}
f{4} = @(x)[ 0, 0, 0, 0,   0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0,   0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0,   0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0,   0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0,   0, -2/5, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 1/4,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;zeros(9,15)];
 

%%%%%%%

%f{5} 

f{5} = @(x) [ (cos(x(1))*sin(x(2)))/cos(x(2)),           -sin(x(1))/(sin(x(2))^2 - 1), 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                  -sin(x(1)),                                  0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
           cos(x(1))/cos(x(2)), -(sin(x(1))*sin(x(2)))/(sin(x(2))^2 - 1), 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                         0,                                  0, 0,   0, 0, 1/6, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                         0,                                  0, 0,   0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                         0,                                  0, 0, 1/4, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0;zeros(9,15)];


%%%%%%%%%

%f{6}


f{6} = @(x) [ -(sin(x(1))*sin(x(2)))/cos(x(2)),           cos(x(1))/cos(x(2))^2, 0,    0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                   -cos(x(1)),                           0, 0,    0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
           -sin(x(1))/cos(x(2)), (cos(x(1))*sin(x(2)))/cos(x(2))^2, 0,    0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                          0,                           0, 0,    0, 1/6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                          0,                           0, 0, -2/5,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                          0,                           0, 0,    0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;zeros(9,15)];



%%%%%%%%

f{7} = @(x) zeros(15,15);
f{8} = @(x) zeros(15,15);
f{9} = @(x) zeros(15,15);
f{10} = @(x) zeros(15,15);
f{11} = @(x) zeros(15,15);
f{12} = @(x) zeros(15,15);
f{13} = @(x) zeros(15,15);
f{14} = @(x) zeros(15,15);
f{15} = @(x) zeros(15,15);

end