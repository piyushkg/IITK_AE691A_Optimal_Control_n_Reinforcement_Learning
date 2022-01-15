
G = 6.673e-11;
AU = 1.5e11;
omega = 1.991e-7;
Msun = 2.0e30;
re = 1.5e11;
rc = 4.5e11;

del_v1 = abs(sqrt((G*Msun)/re)*(sqrt(rc/(rc+re)) - 1))
del_v2 = abs(sqrt((G*Msun)/rc)*(1 - sqrt(re/(rc+re))))
dv = del_v1 + del_v2