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

C = eye(n);

D = zeros(n,1);

sys = ss(A,B,C,D);

%% LQR

Q = diag([7 40 1 10]); R = 0.003;
[K,S,P] = lqr(sys,Q,R,0);

% Defining the closed-loop system
% A_cl = A - B*K;
% B_cl = [0;0;0;0];
% C_cl = C;
% D_cl = D;
% sys_cl = ss(A_cl,B_cl,C_cl, D_cl,0);
% %Setting the initial conditions
% x0 = [0 5*pi/180 0 0]; %Initial rod angle = 5 degrees
% %Plotting the response of the closed-loop system to the initial conditions
% initial(sys_cl,x0)
% grid
% 
% %% Method1
% N = [A B; C D]\[zeros(4); ones(4)];
% Nx = N(1:4,:);Nu = N(5,:);
% 
