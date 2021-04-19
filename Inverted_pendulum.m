clear;

p = parameters;

A = [[0 0 1 0];
    [0 0 0 1];
    [0 -p.m*p.g/p.M -p.K_g^2*p.K_m*p.K_b/(p.M*p.R_m*p.r^2) 0];
    [0 (p.M+p.m)*p.g/(p.m*p.l) p.K_g^2*p.K_m*p.K_b/(p.M*p.R_m*p.r^2*p.l) 0]];
B = [0;
    0;
    p.K_m*p.K_g/(p.M*p.R_m*p.r);
    -p.K_m*p.K_g/(p.r*p.R_m*p.M*p.l)];
C = [[1 0 0 0];
    [0 1 0 0]];
D = 0;