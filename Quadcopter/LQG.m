clear;

load('references_09.mat')
%% Linearization
pq = parameters;
n = 12;%number of states 
m = 4; %number of inputs
p = 6; %number of outputs
pi = 3; %number of integral action outputs
A = zeros(n,n);
A(1,4) = 1;A(2,5) = 1;A(3,6) = 1;
A(4,4) = -pq.kd/pq.m; A(5,5) = -pq.kd/pq.m; A(6,6) = -pq.kd/pq.m;
A(5,7) = -pq.g; A(4,8) = pq.g;
A(7,10) = 1; A(8,11) = 1; A(9,12) = 1;

B = zeros(n,m);
B(6,1) = pq.k*pq.cm/pq.m; B(6,2) = pq.k*pq.cm/pq.m; B(6,3) = pq.k*pq.cm/pq.m; B(6,4) = pq.k*pq.cm/pq.m;
B(10,1) = pq.L*pq.k*pq.cm/pq.Ixx; B(10,3) = -pq.L*pq.k*pq.cm/pq.Ixx;
B(11,2) = pq.L*pq.k*pq.cm/pq.Iyy; B(11,4) = -pq.L*pq.k*pq.cm/pq.Iyy;
B(12,1) = pq.b*pq.cm/pq.Izz; B(12,2) = -pq.b*pq.cm/pq.Izz; B(12,3) = pq.b*pq.cm/pq.Izz; B(12,4) = -pq.b*pq.cm/pq.Izz;

C = zeros(p,n);
C(1:3,1:3) = eye(3); C(4:6,7:9) = eye(3);

D = zeros(p,m);

lin_sys = ss(A,B,C,D);

%% Discretisation

Ts = 0.05;
sysd = c2d(lin_sys, Ts,'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sysd);

%% Integral Action
Ai = [eye(pi) Cd(1:pi,:); zeros(n,pi) Ad];
Bi = [Dd(1:pi,:) ; Bd ];

%% LQR with integral action
% [x y z vx vy vz phi theta psi wx wy wz]
% Q = diag([1 1 3 21 21 1681 6 6 56 110 110 560 16 16 16]); R = 3*diag([1 1 1 1]);
Q = diag([0.7 0.7 0.7 180 180 250 38 38 400 900 900 1400 560 560 700]); 
R = 3e-5*diag([1 1 1 1]);
[K,S,P] = dlqr(Ai,Bi,Q,R,0);

Ki = K(:,1:3);
Ks = K(:,3+1:end);

%% Kalman filter for state estimation

B1 = eye(n);
W = eye(n,n);
% W(1,1) = 2.5*10^-5; W(2,2) = 2.5*10^-5; W(3,3) = 2.5*10^-5;
% W(7,7) = 7.57*10^-5; W(8,8) = 7.57*10^-5; W(9,9) = 7.57*10^-5;
V = [2.5*10^-5*eye(3), zeros(3,3); zeros(3,3), 7.57*10^-5*eye(3)];
[L,P,Z,E] = dlqe(Ad,B1,Cd,W,V);
  

