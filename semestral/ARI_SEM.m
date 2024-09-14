%-----------------------------------------------------
% Linearizace - Symbolické matice A,B,C,D

clc
clear



syms v_T S_T kc h1 h2 g v_O1 v_O S_u S_O S u h_1off h_2off

Inputs.u = u;
Inputs.v_T = v_T;
Inputs.v_O = v_O;

States.h1 = h1;
States.h2 = h2;

Outputs.y1 = h1;
Outputs.y2 = h2;


Symbolic.coupling = Inputs.v_T*S_T*sign(States.h1-States.h2)*sqrt(2*g*abs(States.h1-States.h2));
Symbolic.dh1 = 1/S*(kc*Inputs.u ...
        - v_O1*S_u*sqrt(2*g*(States.h1+h_1off)) ...
        -Symbolic.coupling );
Symbolic.dh2 = 1/S*( Symbolic.coupling ...
        - Inputs.v_O*S_O*sqrt(2*g*(States.h2+h_2off)));

coupling = v_T*S_T*sign(h1-h2)*sqrt(2*g*abs(h1-h2));
dh1 = 1/S*( kc*u - v_O1*S_u*sqrt(2*g*(h1+h_1off)) -coupling );
dh2 = 1/S*( coupling - v_O*S_O*sqrt(2*g*(h2+h_2off)));

disp(dh1,dh2)
Symbolic.A = [diff(Symbolic.dh1, States.h1) diff(Symbolic.dh1, States.h2); diff(Symbolic.dh2, States.h1) diff(Symbolic.dh2, States.h2)];

Symbolic.B = [diff(Symbolic.dh1, Inputs.u), diff(Symbolic.dh1, Inputs.v_T), diff(Symbolic.dh1, Inputs.v_O); diff(Symbolic.dh2, Inputs.u), diff(Symbolic.dh2, Inputs.v_T), diff(Symbolic.dh2, Inputs.v_O)];

Symbolic.C = [diff(Outputs.y1, States.h1) diff(Outputs.y1, States.h2); diff(Outputs.y2, States.h1) diff(Outputs.y2, States.h2)];

Symbolic.D = [diff(Outputs.y1, Inputs.u), diff(Outputs.y1, Inputs.v_T), diff(Outputs.y1, Inputs.v_O);diff(Outputs.y2, Inputs.u), diff(Outputs.y2, Inputs.v_T), diff(Outputs.y2, Inputs.v_O)];


disp('Symbolické matice vypočteny')

clear v_T S_T kc h1 h2 g v_O1 v_O S_u S_O S u h_1off h_2off

%%
% Pracovní bod P

P.u = 2.5;
P.h1 = 0.05182;
P.h2 = 0.01872;
P.v_O = 0.3;
P.v_T = 0.4;

disp('Pracovní bod nastaven')


%%
%-----------------------------------------------------
% Linearizace - Numerické matice An,Bn,Cn,Dn


%p = 1000;
%u_N = 0;


Constants.v_O1 = 0;



S = 50.00/10000;
kc = 3.35E-5;
g = 9.81;
v_O1 = 0;
S_T = 2.6/10000;
S_u = 4.6/10000;
S_O = 4.6/10000;
h_1off = 0;
h_2off = 0;

u = P.u;
h1 = P.h1;
h2 = P.h2;
v_T = P.v_T;
v_O = P.v_O;


An = vpa(subs(Symbolic.A));
Bn = vpa(subs(Symbolic.B));
Cn = vpa(subs(Symbolic.C));
Dn = vpa(subs(Symbolic.D));

P.A = double(An);
P.B = double(Bn);
P.C = double(Cn);
P.D = double(Dn);

disp('Matice vypočteny pro pracovní bod P')

clear v_T S_T kc h1 h2 g v_O1 v_O S_u S_O S u h_1off h_2off An Bn Cn Dn



%%
%-----------------------------------------------------
% RLTool Regulátor - Načtení regulátoru

load C_RL_PI.mat
load SISO_P.mat

disp('Regulátory načteny ze souboru')

%%
%-----------------------------------------------------
% RLTool Regulátor - Nový regulátor

system = ss(P.A, P.B, P.C, P.D);
rltool(system(2,1))

disp('RLTool spuštěn pro nastavení regulátoru')
clear system

%%
%-----------------------------------------------------
% SISOTool Regulátor - Nový regulátor

system = ss(P.A, P.B, P.C, P.D);
sisotool(system(2,1))

disp('SISOTool spuštěn pro nastavení regulátoru')
clear system

%%
%-----------------------------------------------------
% Stavová zpětná vazba - 


s =tf('s');

G_s = P.C*(s*eye(2)-P.A)^-1*P.B+P.D

%step(G_s)
zeros = tzero(G_s)



A2 = P.A;

B2 = [P.B(1,1); 
      P.B(2,1)];

C2 = [P.C(2,1) P.C(2,2)];

D2 = P.D;



Acl = [A2 [0;0];
       -C2 0]

Bcl = [B2;
       0]

Ccl = [C2 0]

Dcl = D2;


T_s = 10;
OS = 0.2;

zeta = -log(OS)/(sqrt(pi^2+log(OS)^2));
omega_n = -log(0.02 * sqrt(1-zeta^2))/(T_s * zeta);
H_s = omega_n^2/(s^2+2*zeta*omega_n*s+omega_n^2);

zpk_H = zpk(H_s)

xd= [zpk_H.P{1}; zeros(1)]

K = place(Acl,Bcl, xd)


clear A2 Acl B2 Bcl C2 Ccl D2 Dcl H_s omega_n OS s T_s xd zeros zeta zpk_H



%%
%FeedForward
s = tf('s');
transfers = tf(system(2,1));


xd = transfers(1)
G_s = P.C*(s*eye(2)-P.A)^-1*P.B+P.D

xd2 = G_S/xd


%%
%Identifikace

syms v_T S_T h1 h2 g S kc u v_O1 v_O S_u h_1off S_O h_2off

coupling = v_T*S_T*sign(h1-h2)*sqrt(2*g*abs(h1-h2));
dh1 = 1/S*(kc*u - v_O1*S_u*sqrt(2*g*(h1+h_1off)) - coupling );
dh2 = 1/S*( coupling - v_O*S_O*sqrt(2*g*(h2+h_2off)));

coupling = v_T*S_T*sign(h1-h2)*sqrt(2*g*abs(h1-h2));
dh1 = 1/S*(kc*u - coupling );
dh2 = 1/S*( coupling - v_O*S_O*sqrt(2*g*(h2+h_2off)));

v_T = 0;
v_O = 0.3;
u = 5;

h1 = 0;
h2 = 0;


dh1 == 0;
dh2 == 0;


g=9.81;
v_O1 = 0;

new = 1/S * ( kc*u - v_O*S_O*sqrt(2*g*(h2+h_2off))) ==0

subs(new)

%S_T = 0.05199997820083527162599894586552*S;
%S_O = (0.03553119189685374257429335376127*S)/(9.81*h_2off + 0.1492101)^(1/2);
%kc = S*6.7*10^-3;

%u = 1/(6.7*10^-3);

S = 50.00/10000;
%kc = 3.35E-5;
g = 9.81;
v_O1 = 0;
%S_T = 2.6/10000;
S_u = 4.6/10000;
%S_O = 4.6/10000;
h_1off = 0;
h_2off = 0;

disp("kc:")
vpa(subs(kc))

disp("S_O:")
vpa(subs(S_O))

disp("S_T:")
vpa(subs(S_T))


%sol = solve(vpa(subs(xd)), kc)
%vpa(subs(dh2))

clear v_T S_T h1 h2 g S kc u v_O1 v_O S_u h_1off S_O h_2off



%%

syms S_T S_u S_O S kc h_1off h_2off u_N

g = 9.81;
p = 1000;

v_T = 0.3;
v_O = 0.3;
u = 5;

h1 = 2.073*10^-1;
h2 = 7.509*10^-2;

dh1 = sign(kc*(u - u_N)^2 - p*g*(h1+h_1off))*S_u/S*(sqrt(2/p*abs(kc*(u - u_N)^2 - p*g*(h1+h_1off)))-v_T*S_T*sqrt(2*g*abs(h1-h2)));        
dh2 = sign(h1-h2)*v_T*S_T/S*(sqrt(2*g*abs(h1-h2)))-v_O*S_O/S*(sqrt(2*g*(h2+h_2off)));


dh1 == 0;
dh2 == 0;

%vpa(subs(dh1))

%vpa(subs(dh1))
%vpa(subs(dh2))

kc = S*6.7*10^-3;
S_T = 0.05199997820083527162599894586552*S;

assume(S > 0);
assume(S_T > 0);
assume(S_O > 0);
assume(kc > 0);

dh1_2 = (S_u*sign(9810.0*h_1off - 1.0*kc*(u_N - 2.5)^2 + 680.0292)*(0.32211551965094760152652497708914*S_T - 1.0*(0.002*abs(9810.0*h_1off - 1.0*kc*(u_N - 2.5)^2 + 680.0292))^(1/2)))/S == 0;
dh2_2 = (0.32211551965094760152652497708914*S_T)/S - (0.4*S_O*(19.62*h_2off + 0.2071872)^(1/2))/S ==0;

dh1_3 = (S_u*sign(9810.0*h_1off - 1.0*kc*(u_N - 2.5)^2 + 760.5693)*(0.32211551965094760152652497708914*S_T - 1.0*(0.002*abs(9810.0*h_1off - 1.0*kc*(u_N - 2.5)^2 + 760.5693))^(1/2)))/S == 0;
dh2_3 = (0.32211551965094760152652497708914*S_T)/S - (0.3*S_O*(19.62*h_2off + 0.3682674)^(1/2))/S == 0;

dh1_4 = (S_u*sign(9810.0*h_1off - 1.0*kc*(u_N - 1.0)^2 + 81.344520000000007445284389007156)*(0.096631001236663177644814481936919*S_T - 1.0*(0.002*abs(9810.0*h_1off - 1.0*kc*(u_N - 1.0)^2 + 81.344520000000007445284389007156))^(1/2)))/S == 0;
dh2_4 = (0.096631001236663177644814481936919*S_T)/S - (0.3*S_O*(19.62*h_2off + 0.058938480000000003132878956257912)^(1/2))/S == 0;

dh1_5 = (S_u*sign(9810.0*h_1off - 1.0*kc*(u_N - 5.0)^2 + 2033.613)*(0.48317327947642136898309672687901*S_T - 1.0*(0.002*abs(9810.0*h_1off - 1.0*kc*(u_N - 5.0)^2 + 2033.613))^(1/2)))/S == 0;
dh2_5 = (0.48317327947642136898309672687901*S_T)/S - (0.3*S_O*(19.62*h_2off + 1.4732658)^(1/2))/S == 0;

%sol = solve([vpa(subs(dh1)),vpa(subs(dh2))])

sol = solve([dh1_2, dh2_2, dh1_3, dh2_3, dh1_4, dh2_4, dh1_5, dh2_5])

%pretty(sol)
%sol
%sol.kc


clear v_T S_T h1 h2 g S kc u v_O S_u h_1off S_O h_2off u_N p
%%

syms S_T S_u S_O S kc h_1off h_2off u_N

v_T = 0;
v_O = 0.3;
u = 2.5;

h1 = 8.375*10^-2;
h2 = 1.256*10^-14;


g=9.81;


coupling = v_T*S_T*sign(h1-h2)*sqrt(2*g*abs(h1-h2));
dh1 = ( kc*u - coupling ) /S;
dh2 = ( coupling - v_O*S_O*sqrt(2*g*(h2+h_2off))) /S;


h_2off = 0;
S_O = 25O/23 * S

dh1 == 1.675 * 10^-2;
dh2 == 2.220 * 10^-15


%S_T = 0.052*S;
%S_O = (0.03553119189685374257429335376127*S)/(9.81*h_2off + 0.1492101)^(1/2);
%kc = S*6.7*10^-3;

%kc = 0.01675*S;



disp("kc:")
vpa(subs(kc))

disp("S_O:")
vpa(subs(S_O))

disp("S_T:")
vpa(subs(S_T))

disp("dh1:")
vpa(subs(dh1))

disp("dh2:")
vpa(subs(dh2))


clear v_T S_T h1 h2 g S kc u v_O S_u h_1off S_O h_2off u_N p dh1 dh2 xd



%%
% Výpočet h_2off
% nastavím stabilní bod a vezmu dh2, za S_O dosadím (23/13)S_T, z dh2
% vytknu S_T/S (...) vnitřek závorky je 0 <- S_T>0

syms h_2off S_T

v_T = 0.4;
v_O = 0.3;
u = 2.5;

h1 = 5.173*10^-2;
h2 = 1.873*10^-2;

g=9.81;

solv_h_2off = (v_T*sign(h1-h2) * sqrt(2*g*abs(h1-h2)) - v_O *(23/13)*sqrt(2*g*(h2+h_2off)) ) == 0;

disp("h_2off:")
vpa(solve(solv_h_2off,h_2off))
 

clear v_T S_T h1 h2 g S kc u v_O S_u h_1off S_O h_2off u_N p solv_h_2off

%%
% Výpočet S_O
% vezmu dh1 + dh2
% najdu rovnovazny stav
% c se odecte
% 



syms S_T S_u S_O S kc h_1off h_2off u_N C

v_O = 0.4;
u = 2.5;

h2 = 1.056*10^-2;


g=9.81;

kc = 6.7*10^-3;

%S_O/S == C

new = kc*u == v_O*C*sqrt(2*g*(h2+h_2off));

vpa(new)

h_2off = 0
vpa(solve(subs(new), C))*50


clear v_T S_T h1 h2 g S kc u v_O S_u h_1off S_O h_2off u_N p new dh2

%%
S = 50/10000;

S_O = (23/250) * S
S_T = (13/250) * S

%%

LQR
Kalman


%%


disp(Ident_Data)

new_data = [Ident_Data(1) Ident_Data(3)];

disp(new_data)

out_data = Ident_Data(1:end,3)-P.h2
in_data = Ident_Data(1:end,1);

%%

system = ss(P.A, P.B, P.C, P.D)
tf(system(2,1))
clear system



%%

[xdd_a, xdd_b, xdd_c, xdd_d] = tf2ss([0.002848], [1, 1.564, 0.1533])

