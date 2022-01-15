inp0 = [1,2,1,2,1,2,5];
Sol = fsolve(@(inp) q4_fun(inp),inp0);

cp1 = Sol(1);
cp2 = Sol(2);

cq1 = Sol(3);
cq2 = Sol(4);

cr1 = Sol(5);
cr2 = Sol(6);

tf = Sol(7);

dt=0.01;tf;t=[0:dt:tf]';m=length(t);
J1= 0.012;
J2=0.026;
J3=0.038;

 
theta_q = zeros(m,1);
q = zeros(m,1);

for i=1:m,
    theta_q(i,1) = (cq1/J2)*((t(i,1)^3)/6) - (cq2/J2)*((t(i,1)^2)/2) + 180;
    q(i,1) = (cq1/J2)*((t(i,1)^2)/2) - (cq2/J2)*t(i,1);
end;
figure(1)
plot(theta_q,q, 'linewidth', 1.8)
title('Assgn - 2 Q4: Pitch Rate vs Pitch Angle')
xlabel('theta')
ylabel('rate')

figure(2)
plot(t,theta_q, 'linewidth', 1.8)
title('Assgn - 2 Q4: Pitch Angle vs Time')
xlabel('time')
ylabel('theta')

figure(3)
plot(t,q, 'linewidth', 1.8)
title('Assgn - 2 Q4: Pitch rate vs Time')
xlabel('time')
ylabel('rate')

ct1 = zeros(m,1);
ct2 = zeros(m,1);
ct3 = zeros(m,1);
ct4 = zeros(m,1);
m_torq = zeros(m,1);
inp20 = [1,2,1,2];


for i=1:m,
    m_torq(i,1) = (cq1/J2)*t(i,1) - (cq2/J2);
    Sol2 = fsolve(@(inp2) q4_fun2(inp2, m_torq(i,1), theta_q(i,1)),inp20);
    ct1(i,1) = Sol2(1);
    ct2(i,1) = Sol2(2);
    ct3(i,1) = Sol2(3);
    ct4(i,1) = Sol2(4);
end;


figure(4)


plot(t,ct1,'DisplayName', 'ct1','linewidth', 1.8)
hold on

plot(t,ct3, 'DisplayName', 'ct3','linewidth', 1.8)
ylim([-0.1 0.1])
title('Variation of CT1 and CT3 with time')
xlabel('Time')
ylabel('CT Value')
hold off

figure(5)
plot(t,ct2, 'DisplayName', 'ct2','linewidth', 1.8)
hold on
plot(t,ct4, 'DisplayName', 'ct4','linewidth', 1.8)
ylim([-0.1 0.1])
title('Variation of CT2 and CT4 with time')
xlabel('Time')
ylabel('CT Value')
hold off


