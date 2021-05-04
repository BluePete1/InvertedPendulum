% Run first lin_dis.m

%Poles
poles = eig(Ad);
pzmap(Z)

%Controlability
CO = ctrb(Ad,Bd);
rank(CO)
%->Controllable

%PBH
[V, d] = eig(Ad');
Pc = Bd'*V;
%->Stabilizable

%Observable
OB = obsv(Ad,Cd);
rank(OB)
%->Observable

%PBH
[V,d]=eig(Ad);
Po = Cd*V;
%->Detectable

%->Minimal

%transmission zeros
tz = tzero(Ad,Bd,Cd,Dd);
