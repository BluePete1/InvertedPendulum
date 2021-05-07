clear;

p = parameters;
%% Linear state space model

n=4;

A = [[0 0 1 0];
    [0 0 0 1];
    [0 -p.m*p.g/p.M -p.Kg^2*p.Km*p.Kb/(p.M*p.Rm*p.r^2) 0];
    [0 (p.M+p.m)*p.g/(p.M*p.l) p.Kg^2*p.Km*p.Kb/(p.M*p.Rm*p.r^2*p.l) 0]];

B = [0;
    0;
    p.Km*p.Kg/(p.M*p.Rm*p.r);
    -p.Km*p.Kg/(p.r*p.Rm*p.M*p.l)];

C = zeros(n,n); C(1,1) = 1; C(2,2) = 1;

D = zeros(n,1);

sys = ss(A,B,C,D);

%% Measurment inacuracies

wc = 2*2*pi; % 2 Hz
Ts = 1/200;  % 5 ms

V_max_in = 5;
V_max_out = 10;

% Conversion to Volts:
% x_measured = (4.41/0.456)*x
% alfa_measured = (6.328/(pi/4))*alfa
Con_to_V = [4.41/0.456, 0, 0, 0; 0, 6.328/(pi/4), 0, 0];
Con_from_V = [0.456/4.41, 0; 0, (pi/4)/6.328];

% Quantization step
quant_step = 20/(2^16);

%% LQR

% Q = diag([0.25 4 0.0 0.0]); R = 0.003;
Q = diag([7 40 1.5 10]); R = 0.003;
% Q = diag([1 4 0.0 0.0]); R = 0.003;
[K,S,P] = lqr(sys,Q,R,0);

% Defining the closed-loop system
A_cl = A - B*K;
B_cl = B*K; % Nieuwe input is x_ref
C_cl = C-D*K;
D_cl = D*K; % Nieuwe input is x_ref
sys_cl = ss(A_cl,B_cl,C_cl, D_cl);

%Setting the initial conditions
x0 = [0 5*pi/180 0 0]; %Initial rod angle = 5 degrees

% Closed-loop eigenvalues
lambdas = eig(A_cl);

figure
plot(real(lambdas),imag(lambdas),'bx')
grid
xlabel('Re')
ylabel('Im')

% Step response
figure
step(sys_cl(:,1))

% Impulse response
figure
impulse(sys_cl(:,1))

% Begin at angle
figure
initial(sys_cl,x0)
grid
















