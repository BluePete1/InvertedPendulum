clear;

load('references_09.mat')
%% Linearization
pq = parameters;
n = 12;%number of states 
m = 4; %number of inputs
p = 6; %number of outputs
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
C(1,1) = 1; C(2,2) = 1; C(3,3) = 1;
C(4,7) = 1; C(5,8) = 1; C(6,9) = 1;
% C = eye(p,n);

D = zeros(p,m);

lin_sys = ss(A,B,C,D);

%% Discretisation

Ts = 0.05;
sysd = c2d(lin_sys, Ts,'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sysd);

%% Full state feedback controller
big_A = [Ad-eye(size(Ad)) Bd
         eye(12,12) zeros(12,m)];

big_Y =[ zeros(n,12)
         eye(12,12) ];

big_N = (big_A)\big_Y;

Nx = big_N(1:n,:);
Nu = big_N (n+1:end,:);

%% Pole placement Contoller

P = [0.94 0.94 0.94 0.67 0.67 0.67 0.92 0.92 0.92 0.82 0.82 0.82];

K = place(Ad,Bd,P);
%% Pole placement Estimator

L = place(Ad',Cd',P)';
Dpp = zeros(12,10);

