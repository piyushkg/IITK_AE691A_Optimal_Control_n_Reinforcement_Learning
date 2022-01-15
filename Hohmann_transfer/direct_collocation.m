
import casadi.*
clc;
% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = [0 collocation_points(d, 'legendre')];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:d+1
    if r ~= j
      coeff = conv(coeff, [1, -tau_root(r)]);
      coeff = coeff / (tau_root(j)-tau_root(r));
    end
  end
  % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D(j) = polyval(coeff, 1.0);

  % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);
  for r=1:d+1
    C(j,r) = polyval(pder, tau_root(r));
  end

  % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);
  B(j) = polyval(pint, 1.0);
end

% Time horizon
T = 473;

% Declare model variables
x = SX.sym('x');
y = SX.sym('y');
vx = SX.sym('vx');
vy = SX.sym('vy');
state = [x;y;vx;vy];
ax= SX.sym('ax');
ay = SX.sym('ay');
r = sqrt(x^2 + y^2);
mu_con = 2.95057084 *10^-4;
u = [ax;ay];

% Model equations
% state_dot = [vx; vy; ((-mu_con*x)/ (r^3)) + ax;  ((-mu_con*y) / (r^3)) + ay];
state_dot = [vx; vy; ((-mu_con*x)/ ((sqrt(x^2 + y^2))^3)) + ax;  ((-mu_con*y) / ((sqrt(x^2 + y^2))^3)) + ay];

% Objective term
L = 0.5*(ax^2 + ay^2);

% Continuous time dynamics
f = Function('f', {state, u}, {state_dot, L});

% Control discretization
N = 200; % number of control intervals
h = T/N;

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
Xk = MX.sym('X0', 4);
w = {w{:}, Xk};
lbw = [lbw; 0; -1;0.01728; 0];
ubw = [ubw; 0; -1;0.01728; 0];
w0 = [w0; 0; -1;0.01728; 0];

% Formulate the NLP
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],2);
    w = {w{:}, Uk};
    lbw = [lbw; -inf; -inf];
    ubw = [ubw; inf; inf];
    w0 = [w0; 0; 0];

    % State at collocation points
    Xkj = {};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], 4);
        w = {w{:}, Xkj{j}};
        lbw = [lbw; -inf;-inf;-inf;-inf];
        ubw = [ubw;  inf;inf;inf;inf];
        w0 = [w0; 0; -1;0.01728; 0];
    end

    % Loop over collocation points
    Xk_end = D(1)*Xk;
    for j=1:d
       % Expression for the state derivative at the collocation point
       xp = C(1,j+1)*Xk;
       for r=1:d
           xp = xp + C(r+1,j+1)*Xkj{r};
       end

       % Append collocation equations
       [fj, qj] = f(Xkj{j},Uk);
       g = {g{:}, h*fj - xp};
       lbg = [lbg; 0; 0;0;0];
       ubg = [ubg; 0; 0; 0;0];

       % Add contribution to the end state
       Xk_end = Xk_end + D(j+1)*Xkj{j};

       % Add contribution to quadrature function
       J = J + B(j+1)*qj*h;
    end

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], 4);
    w = {w{:}, Xk};
    if(k==N-1)
        lbw = [lbw;0; 2.766;-0.10368;0.0];
        ubw = [ubw;0; 2.766;-0.10368;0.0];
        w0 = [w0; 0; 2.766;-0.010368; 0];
    else
        lbw = [lbw;-inf;-inf;-inf;-inf];
        ubw = [ubw;  inf;inf;inf;inf];
        w0 = [w0; 0; 0;0; 0];
    end

    % Add equality constraint
    g = {g{:}, Xk_end-Xk};
    lbg = [lbg;  0; 0;0;0];
    ubg = [ubg;  0; 0;0;0];
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);
size(w_opt)

% Plot the solution
x1_opt = w_opt(1:6+4*d:end);
x2_opt = w_opt(2:6+4*d:end);
x3_opt = w_opt(3:6+4*d:end);
x4_opt = w_opt(4:6+4*d:end);
u1_opt = w_opt(5:6+4*d:end);
u2_opt = w_opt(6:6+4*d:end);
tgrid = linspace(0, T, N+1);
clf;
hold on
size(tgrid)
size(x1_opt)
plot(tgrid, x1_opt', '--')
plot(tgrid, x2_opt', '-')
xlabel('t')
% legend('x1','x2','u')
