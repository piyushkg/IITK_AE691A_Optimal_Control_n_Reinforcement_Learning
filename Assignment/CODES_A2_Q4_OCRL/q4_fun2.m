function fval = q4_fun2(inp, m_inst, theta)
ct1 = inp(1);
ct2 = inp(2);
ct3 = inp(3);
ct4 = inp(4);

nu = 1;
m = 1.4;
R = 0.14;
g = 9.8;
l = 0.18;
V = R*418.9;
ro=1;
k = ro*pi*(R^2)*(V^2);
cdtheta = cosd(theta);

% if cdtheta < 0.5
%     cdtheta = 0.5;
% elseif cdtheta > -0.5
%     cdtheta = -0.5;
    
    
fval(1,1) =  nu*k*(ct1 + ct2 + ct3 + ct4) - (m*g/cdtheta) ;
fval(2,1) = nu*k*l*(ct1 - ct2 - ct3 + ct4);

fval(3,1) =  nu*k*l*(ct1 + ct2 - ct3 - ct4) - m_inst ;
fval(4,1) = nu*k*R*(1/sqrt(2))*((abs(ct1)^1.5) - (abs(ct2)^1.5) + (abs(ct3)^1.5) - (abs(ct4)^1.5));

end