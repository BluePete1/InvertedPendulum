clear;

load('references_09.mat')
%% Linearization
pq = parameters;
n = 12;%number of states 
m = 4; %number of inputs
p = 12; %number of outputs
A = zeros(n,n);
A(1,4) = 1;A(2,5) = 1;A(3,6) = 1;
A(4,4) = -pq.kd/pq.m; A(5,5) = -pq.kd/pq.m; A(6,6) = -pq.kd/pq.m;
A(5,7) = -pq.g; A(4,8) = pq.g;
A(7,10) = 1; A(8,11) = 1; A(9,12) = 1;

B = zeros(12,4);
B(6,1) = pq.k*pq.cm/pq.m; B(6,2) = pq.k*pq.cm/pq.m; B(6,3) = pq.k*pq.cm/pq.m; B(6,4) = pq.k*pq.cm/pq.m;
B(10,1) = pq.L*pq.k*pq.cm/pq.Ixx; B(10,3) = -pq.L*pq.k*pq.cm/pq.Ixx;
B(11,2) = pq.L*pq.k*pq.cm/pq.Iyy; B(11,4) = -pq.L*pq.k*pq.cm/pq.Iyy;
B(12,1) = pq.b*pq.cm/pq.Izz; B(12,2) = -pq.b*pq.cm/pq.Izz; B(12,3) = pq.b*pq.cm/pq.Izz; B(12,4) = -pq.b*pq.cm/pq.Izz;

% C = zeros(12,12);
% C(1,1) = 1; C(2,2) = 1; C(3,3) = 1;
% C(7,7) = 1; C(8,8) = 1; C(9,9) = 1;

C = eye(n);

D = zeros(n,m);

lin_sys = ss(A,B,C,D);

%% Discretisation

Ts = 0.05;
sysd = c2d(lin_sys, Ts,'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sysd);

%% LQR
% [x y z vx vy vz phi theta psi wx wy wz]
Q = diag([21 14 81 6 6 6 110 160 560 9 6 6]); R = diag([1 1 1 1]);
[K,S,P] = dlqr(Ad,Bd,Q,R,0);

big_A = [Ad-eye(size(Ad)) Bd
         Cd Dd];

big_Y =[ zeros(n,p)
         eye(p,p) ];

big_N = (big_A)\big_Y;

Nx = big_N(1:n,:)
Nu = big_N (n+1:end,:)
