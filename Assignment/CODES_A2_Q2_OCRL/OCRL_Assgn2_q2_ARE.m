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
[P,L,G] = care(A,B,Q,R)
X1 = zeros(6,data_size(1)+1);
X1(:,1) = X0;

X1_ARE = zeros(6,data_size(1)+1);
X1_ARE(:,1) = X0;

% Solving close loop dynamical equation
for i = 1:1:data_size(1)
    P_inst =  zeros(6,6)
    P_inst(:,1) = P_flip(i,1:6);
    P_inst(:,2) = P_flip(i,7:12);
    P_inst(:,3) = P_flip(i,13:18);
    P_inst(:,4) = P_flip(i,19:24);
    P_inst(:,5) = P_flip(i,25:30);
    P_inst(:,6) = P_flip(i,31:36);
    X_dot = (A-B*(R_inv*B'*P_inst))*X1(:,i);
    X1(:,i+1) = X1(:,i) + dt*X_dot;
    
    X_dot_ARE = (A-B*(R_inv*B'*P))*X1_ARE(:,i);
    X1_ARE(:,i+1) = X1_ARE(:,i) + dt*X_dot_ARE;
end

%% CASE-2
% Define Q and R Matrices
q1 = 100;q2 = 100;q3 = 1;q4 = 1;q5 = 1;q6 = 1;
r1 = 1;r2 = 100;r3 = 1;
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
[P,L,G] = care(A,B,Q,R)
X2 = zeros(6,data_size(1)+1);
X2(:,1) = X0;

X2_ARE = zeros(6,data_size(1)+1);
X2_ARE(:,1) = X0;

% Solving close loop dynamical equation
for i = 1:1:data_size(1)
    P_inst =  zeros(6,6)
    P_inst(:,1) = P_flip(i,1:6);
    P_inst(:,2) = P_flip(i,7:12);
    P_inst(:,3) = P_flip(i,13:18);
    P_inst(:,4) = P_flip(i,19:24);
    P_inst(:,5) = P_flip(i,25:30);
    P_inst(:,6) = P_flip(i,31:36);
    X_dot = (A-B*(R_inv*B'*P_inst))*X2(:,i);
    X2(:,i+1) = X2(:,i) + dt*X_dot;
    
    X_dot_ARE = (A-B*(R_inv*B'*P))*X2_ARE(:,i);
    X2_ARE(:,i+1) = X2_ARE(:,i) + dt*X_dot_ARE;
end

figure(1)
plot(t_for,X1(2,1:10001),'DisplayName', 'Case-1-DRE (S1=1, S2=1)','linewidth', 1.8)
hold on
plot(t_for,X1_ARE(2,1:10001),'DisplayName', 'Case-1-ARE (S1=1, S2=1)','linewidth', 1.8)

title('Roll Rate Convergence for Case1 DRE vs ARE')
xlabel('Time')
ylabel('Roll Rate')
hold off

figure(2)
plot(t_for,X1(2,1:10001),'DisplayName', 'Case-2-DRE (S1=100, S2=100)','linewidth', 1.8)
hold on
plot(t_for,X1_ARE(2,1:10001),'DisplayName', 'Case-2-ARE (S1=100, S2=100)','linewidth', 1.8)

title('Roll Rate Convergence for Case2 DRE vs ARE')
xlabel('Time')
ylabel('Roll Rate')
hold off

