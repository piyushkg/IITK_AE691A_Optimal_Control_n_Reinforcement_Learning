% Define various matrices 
A = [0 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 0 0];
B = [0 0 0; 1 0 0 ; 0 0 0; 0 1 0; 0 0 0; 0 0 1];

diag_Sf = [1/s1^2, 1/s2^2, 1/s3^2, 1/s4^2, 1/s5^2, 1/s6^2];
Sf = diag(diag_Sf);
Sf = reshape(Sf,[], 1);

% Define Time Parameters
t0 = 0;
dt = 0.001
tf = 10;
t_for = [t0:dt:tf]';
t_rev = flip(t_for);

% Initial State
X0 = [50; 0; 10; 0; 25; 0];

%% CASE-1
% Define Q and R Matrices
q1 = 1;q2 = 1;q3 = 1;q4 = 1;q5 = 1;q6 = 1;
r1 = 1;r2 = 1;r3 = 1;
s1 = 1;s2 = 1;s3 = 1;s4 = 1;s5 = 1;s6 = 1;

diag_Q = [1/q1^2, 1/q2^2, 1/q3^2, 1/q4^2, 1/q5^2, 1/q6^2];
Q = diag(diag_Q);

diag_R = [1/r1^2, 1/r2^2, 1/r3^2];
R = diag(diag_R);
R_inv = inv(R)

% Solving Differential Ricatti Equation
[t_sol, P_mat] = ode45(@(t,P)P_fun(t,P, A, B, Q, R), t_rev, Sf);
P_flip = flip(P_mat);
data_size = size(P_flip);
X1 = zeros(6,data_size(1)+1);
X1(:,1) = X0;

% Solving close loop dynamical equation
for i = 1:1:data_size(1)
    P_inst = zeros(6,6)
    P_inst(:,1) = P_flip(i,1:6);
    P_inst(:,2) = P_flip(i,7:12);
    P_inst(:,3) = P_flip(i,13:18);
    P_inst(:,4) = P_flip(i,19:24);
    P_inst(:,5) = P_flip(i,25:30);
    P_inst(:,6) = P_flip(i,31:36);
    X_dot = (A-B*(R_inv*B'*P_inst))*X1(:,i);
    X1(:,i+1) = X1(:,i) + dt*X_dot;
end

%% CASE-2
% Define Q and R Matrices
q1 = 100;q2 = 100;q3 = 1;q4 = 1;q5 = 10;q6 = 10;
r1 = 1;r2 = 10;r3 = 1;
s1 = 1;s2 = 1;s3 = 1;s4 = 1;s5 = 1;s6 = 1;

diag_Q = [1/q1^2, 1/q2^2, 1/q3^2, 1/q4^2, 1/q5^2, 1/q6^2];
Q = diag(diag_Q);

diag_R = [1/r1^2, 1/r2^2, 1/r3^2];
R = diag(diag_R);
R_inv = inv(R)

% Solving Differential Ricatti Equation
[t_sol, P_mat] = ode45(@(t,P)P_fun(t,P, A, B, Q, R), t_rev, Sf);
P_flip = flip(P_mat);
data_size = size(P_flip);
X2 = zeros(6,data_size(1)+1);
X2(:,1) = X0;

% Solving close loop dynamical equation
for i = 1:1:data_size(1)
    P_inst = zeros(6,6)
    P_inst(:,1) = P_flip(i,1:6);
    P_inst(:,2) = P_flip(i,7:12);
    P_inst(:,3) = P_flip(i,13:18);
    P_inst(:,4) = P_flip(i,19:24);
    P_inst(:,5) = P_flip(i,25:30);
    P_inst(:,6) = P_flip(i,31:36);
    X_dot = (A-B*(R_inv*B'*P_inst))*X2(:,i);
    X2(:,i+1) = X2(:,i) + dt*X_dot;
end

%% CASE-3
% Define Q and R Matrices
q1 = 1;q2 = 1;q3 = 1;q4 = 1;q5 = 100;q6 = 100;
r1 = 100;r2 = 100;r3 = 1;
s1 = 1;s2 = 1;s3 = 1;s4 = 1;s5 = 1;s6 = 1;

diag_Q = [1/q1^2, 1/q2^2, 1/q3^2, 1/q4^2, 1/q5^2, 1/q6^2];
Q = diag(diag_Q);

diag_R = [1/r1^2, 1/r2^2, 1/r3^2];
R = diag(diag_R);
R_inv = inv(R)

% Solving Differential Ricatti Equation
[t_sol, P_mat] = ode45(@(t,P)P_fun(t,P, A, B, Q, R), t_rev, Sf);
P_flip = flip(P_mat);
data_size = size(P_flip);
X3 = zeros(6,data_size(1)+1);
X3(:,1) = X0;

% Solving close loop dynamical equation
for i = 1:1:data_size(1)
    P_inst = zeros(6,6)
    P_inst(:,1) = P_flip(i,1:6);
    P_inst(:,2) = P_flip(i,7:12);
    P_inst(:,3) = P_flip(i,13:18);
    P_inst(:,4) = P_flip(i,19:24);
    P_inst(:,5) = P_flip(i,25:30);
    P_inst(:,6) = P_flip(i,31:36);
    X_dot = (A-B*(R_inv*B'*P_inst))*X3(:,i);
    X3(:,i+1) = X3(:,i) + dt*X_dot;
end

figure(1)
plot(t_for,X1(1,1:10001),'DisplayName', 'Case-1','linewidth', 1.8)
hold on
plot(t_for,X2(1,1:10001),'DisplayName', 'Case-2','linewidth', 1.8)
hold on
plot(t_for,X3(1,1:10001),'DisplayName', 'Case-3','linewidth', 1.8)
title('Roll Angle Convergence for Case 1-2-3')
xlabel('Time')
ylabel('Roll Angle')
hold off

figure(2)
plot(t_for,X1(2,1:10001),'DisplayName', 'Case-1','linewidth', 1.8)
hold on
plot(t_for,X2(2,1:10001),'DisplayName', 'Case-2','linewidth', 1.8)
hold on
plot(t_for,X3(2,1:10001),'DisplayName', 'Case-3','linewidth', 1.8)
title('Roll Rate Convergence for Case 1-2-3')
xlabel('Time')
ylabel('Roll Rate')
hold off

figure(3)
plot(t_for,X1(3,1:10001),'DisplayName', 'Case-1','linewidth', 1.8)
hold on
plot(t_for,X2(3,1:10001),'DisplayName', 'Case-2','linewidth', 1.8)
hold on
plot(t_for,X3(3,1:10001),'DisplayName', 'Case-3','linewidth', 1.8)
title('Pitch Angle Convergence for Case 1-2-3')
xlabel('Time')
ylabel('Pitch Angle')
hold off

figure(4)
plot(t_for,X1(4,1:10001),'DisplayName', 'Case-1','linewidth', 1.8)
hold on
plot(t_for,X2(4,1:10001),'DisplayName', 'Case-2','linewidth', 1.8)
hold on
plot(t_for,X3(4,1:10001),'DisplayName', 'Case-3','linewidth', 1.8)
title('Pitch Rate Convergence for Case 1-2-3')
xlabel('Time')
ylabel('Pitch Rate')
hold off

figure(5)
plot(t_for,X1(5,1:10001),'DisplayName', 'Case-1','linewidth', 1.8)
hold on
plot(t_for,X2(5,1:10001),'DisplayName', 'Case-2','linewidth', 1.8)
hold on
plot(t_for,X3(5,1:10001),'DisplayName', 'Case-3','linewidth', 1.8)
title('Yaw Angle Convergence for Case 1-2-3')
xlabel('Time')
ylabel('Yaw Angle')
hold off

figure(6)
plot(t_for,X1(6,1:10001),'DisplayName', 'Case-1','linewidth', 1.8)
hold on
plot(t_for,X2(6,1:10001),'DisplayName', 'Case-2','linewidth', 1.8)
hold on
plot(t_for,X3(6,1:10001),'DisplayName', 'Case-3','linewidth', 1.8)
title('Yaw Rate Convergence for Case 1-2-3')
xlabel('Time')
ylabel('Yaw Rate')
hold off