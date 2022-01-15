clear all
close all
%% Defining System Parameters
dtime = 0.1;  % Time step
N = 20;       % Number of time steps
% State Transition Matrix (A_k)
A_k =  [1 0 dtime   0;
        0 1   0   dtime;
        0 0   1     0;
        0 0   0     1];
% Define B_k matrix
B_k = [1*(dtime^2)      0;
        0        1*(dtime^2);
        dtime           0;
        0           dtime];
% Define D_k matrix
eta = 0.01;
D_k = eta*eye(4);  

% Initializing pre-augmented A, B and D matrices
preaugA = A_k; 
preaugB = B_k;
preaugD = D_k;
dim_u = length(B_k(1,:));

% Define Q and R matrices
Q_0 =  [0.5000     0         0         0;
         0    4.0000         0         0;
         0         0    0.0500         0;
         0         0         0    0.0500];
R_0 =  [20  0;
        0  20];
    
augQ = blkdiag(Q_0);
augR  = [];

% Define Initial Conditions
mu_0  = [-10;1;0;0]; 
sig_0 =  9*[0.1    0         0         0;
           0     0.1         0         0;
           0       0      0.01         0;
           0       0         0      0.01];
mu_f  = [0;0;0;0];
sig_f = 0.2*sig_0;      

%% System Initialization and Augmented Matrix formation
global dim_x
global dim_xfull
dim_x = length(mu_0);
dim_xfull = (N+1)*dim_x;

w = randn(dim_x,1); 
X_0 = mu_0 + preaugD*w;

A = eye(dim_x); B = zeros(dim_x,dim_u); D = D_k;
X = X_0; X_bar = mu_0; X_dev_0 = X - X_bar;

for i=1:N
    A = [A; preaugA];
    B = [[B, zeros(length(B(:,1)), length(preaugB(1,:))-length(B(1,:)))]; preaugB];
    D = [[D, zeros(length(D(:,1)), length(preaugD(1,:))-length(D(1,:)))]; preaugD];
    if i<N
        [preaugA,preaugB,preaugD] = create_aug_mat(preaugA,preaugB,preaugD, A_k,B_k,D_k);
    end
    augQ = blkdiag(augQ, Q_0);
    augR = blkdiag(augR, R_0);
end
augA = A;
augB = B;
augD = D;

%% Mean Steering
% Mean Optimization
cvx_begin quiet
    variables V(dim_u*N)    % V - Optimization Variable
    minimize((augA*mu_0+augB*V)'*augQ*(augA*mu_0+augB*V) + V'*augR*V)
    x = augB*V + augA*mu_0;
    subject to
        E_pick(0)*x == mu_0;
        E_pick(N+1)*x == mu_f;
cvx_end
U_bar = V;

%% Close Loop Mean Steering
% R = augB'*augQ*augB+augR;
% U_bar = (R^-1)*(augB'*augQ*augA*mu_0 + preaugB'*((preaugB*(R^-1)*preaugB')^-1)*(mu_f - preaugA*mu_0 - preaugB*(R^-1)*augB'*augQ*augA*mu_0));

%% Covariance Steering
% Some pre-requisites
M = augA*sig_0*augA' + augD*augD';
[V,U]=eig(M); U_sqrt = sqrt(U); M_sqrt = V*U_sqrt*V';
Q_sqrt = sqrt(augQ);
L = [];

% Without Covariance steering section starts
% L = zeros(N*dim_u, N*dim_x); [n_row,n_col] = size(L);
% L = [L,zeros(n_row,dim_x)];
% Without Covariance steering section ends

% Covariance Optimization
cvx_begin 
    variable k(dim_u,dim_x,N);
    for i=1:N
        L = blkdiag(L,k(:,:,i));    % L - Optimization Variable
    end
    [n_row,n_col] = size(L);
    L = [L,zeros(n_row,dim_x)];
    minimize(norm((augQ^0.5)*(eye(dim_x*(N+1))+augB*L)*V*U_sqrt,'fro') + norm((augR^0.5)*L*V*U_sqrt,'fro'))
    subject to
          1-norm(M_sqrt*((eye(dim_x*(N+1))+augB*L)')*(E_pick(N+1)')*(sig_f^(-0.5)))>=0;
cvx_end

%% Post Optimization: State and Control Computation 
% Noise
W = [w;randn((N-1)*dim_x,1)];
% State
X_bar = augA*mu_0 + augB*U_bar;
X_dev = (eye(dim_x*(N+1))+augB*L)*(augA*X_dev_0 + augD*W);
X = X_bar + X_dev;
% Control
U_dev = L*(augA*X_dev_0 + augD*W);
U = U_bar + U_dev;
% Overall Covariance Matrix
cov_full = (eye(dim_x*(N+1))+augB*L)*M*(eye(dim_x*(N+1))+augB*L)';
% Prepare Actual state and Mean state for plotting
X_act = []; Y_act = [];
X_mean = []; Y_mean = [];
for i=1:dim_x:dim_x*(N+1)
    X_mean = [X_mean;X_bar(i)]; Y_mean = [Y_mean;X_bar(i+1)];
    X_act = [X_act;X(i)]; Y_act = [Y_act;X(i+1)];
end

%% Plotting results 
% Plotting Initial and Final Position and Covariance Ellipses 
mu_pos_0 = mu_0(1:2);
cov_pos_0 = sig_0(1:2,1:2);
mu_pos_f = mu_f(1:2);
cov_pos_f = sig_f(1:2,1:2);

figure
plot_cov_elps(mu_pos_0,cov_pos_0,'r', 2)
hold on;
plot_cov_elps(mu_pos_f,cov_pos_f,'r', 2)

%Plotting Actual state and Mean state
plot_x = plot(X_act,Y_act);
plot_xmean = plot(X_mean,Y_mean);
legend([plot_x plot_xmean],'X','X mean')
xlabel('x position'); ylabel('y position'); title('Covariance Steering')

%Plotting Covariance Matrix for each time step
lw=1;
for i=1:dim_x:dim_x*(N+1)
    plot_cov_elps(X_bar(i:i+1),cov_full(i:i+1,i:i+1),'k', lw)
end
hold off





