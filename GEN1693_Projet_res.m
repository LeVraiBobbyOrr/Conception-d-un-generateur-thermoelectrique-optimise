clear;clc
% Constantes
GEN1693_Projet_var
% x = 0.5/1000;

% Fonctionnelles [L,H,x]
fun = @(X) -(((floor((Ltot - X(3))./(X(1) + X(3)))).^2)*alp^2.*((K./(K+2*(2*((floor((Ltot - X(3))./(X(1) + X(3)))).^2)*kTE.*X(1).^2*Kp./(2*((floor((Ltot - X(3))./(X(1) + X(3)))).^2)*kTE.*X(1).^2 + Kp*X(2))))).* (mcp*Theta./((((2*((floor((Ltot - X(3))./(X(1) + X(3)))).^2)*kTE.*X(1).^2*Kp./(2*((floor((Ltot - X(3))./(X(1) + X(3)))).^2)*kTE.*X(1).^2 + Kp*X(2)))*K)/(K+(2*((floor((Ltot - X(3))./(X(1) + X(3)))).^2)*kTE.*X(1).^2*Kp./(2*((floor((Ltot - X(3))./(X(1) + X(3)))).^2)*kTE.*X(1).^2 + Kp*X(2)))))+mcp))).^2./(4*pTE*(X(2)./X(1).^2)))./(X(1).^2*X(2)*(floor((Ltot - X(3))./(X(1) + X(3)))).^2);

% Point de départ
x0 = [1.5/1000 1.5/1000 0.9/1000];

% Options d'optimisation
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','OptimalityTolerance',1e-10);

% Résolution du problème 
[X,fonctionnelle,exitflag,output] = fmincon(fun,x0,[],[],[],[],[1.2/1000, 1/1000, 1/1000],[1.4/1000, 2/1000, 2/1000],[],options);

% Valeur de L,H et x
X

% Valeur fonctionnelle
-fonctionnelle

% Nombre de modules
L = X(1);
H = X(2);
x = X(3);
N = (floor((Ltot - x)./(L + x))).^2;

formatSpec = 'Longueur : %4.2f mm\nHauteur %4.2f mm\nEspace entre module : %4.2f mm\n\nPuissance volumique : %4.2f W/m^2\n';
fprintf(formatSpec, L*1000,H*1000,x*1000,-fonctionnelle);

formatSpec = 'Nombre de module : %4.0f\n';
fprintf(formatSpec, N);
