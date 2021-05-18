clear;

p = parameters;
%% Linear state space model

A = [[0 0 1 0];
    [0 0 0 1];
    [0 -p.m*p.g/p.M -p.Kg^2*p.Km*p.Kb/(p.M*p.Rm*p.r^2) 0];
    [0 (p.M+p.m)*p.g/(p.M*p.l) p.Kg^2*p.Km*p.Kb/(p.M*p.Rm*p.r^2*p.l) 0]];

B = [0;
    0;
    p.Km*p.Kg/(p.M*p.Rm*p.r);
    -p.Km*p.Kg/(p.r*p.Rm*p.M*p.l)];

C = [4.41/0.456, 0, 0, 0;
    0, 6.328/(pi/4), 0, 0];

D = [0;0];

sys = ss(A,B,C,D);

% %% Open loop analysis
% 
% %Poles
% poles = eig(A);
% 
% %transmission zeros
% tz = tzero(A,B,C,D);
% 
% %Controlability
% CO = ctrb(A,B);
% rank(CO)
% %->Controllable
% 
% %PBH
% [V, dc] = eig(A');
% Pc = B'*V;
% %->Stabilizable
% 
% %Observable
% OB = obsv(A,C);
% rank(OB)
% %->Observable
% 
% %PBH
% [V,do]=eig(A);
% Po = C*V;
% %->Detectable
% 
% %->Minimal
