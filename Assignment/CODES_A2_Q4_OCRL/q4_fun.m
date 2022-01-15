function fval = q4_fun(inp)
cp1 = inp(1);
cp2 = inp(2);

cq1 = inp(3);
cq2 = inp(4);

cr1 = inp(5);
cr2 = inp(6);

tf = inp(7);

J1= 0.012;
J2=0.026;
J3=0.038;

fval(1,1) =  (cp1/J1)*((tf^2)/2) - (cp2/J1)*tf ;
fval(2,1) = (cp1/J1)*((tf^3)/6) - (cp2/J1)*((tf^2)/2);

fval(3,1) =  (cq1/J2)*((tf^2)/2) - (cq2/J2)*tf ;
fval(4,1) = (cq1/J2)*((tf^3)/6) - (cq2/J2)*((tf^2)/2) + 180;

fval(5,1) =  (cr1/J3)*((tf^2)/2) - (cr2/J3)*tf ;
fval(6,1) = (cr1/J3)*((tf^3)/6) - (cr2/J3)*((tf^2)/2);

fval(7,1) = 4*(J2^2)*tf - (cq1^2)*(tf^2) + 2*cq1*cq2*tf - (cq2^2)*(tf^2) ;

end