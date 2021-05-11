clear all;
clc;

GEN1693_Projet_var

% Variables : H (hateur element ), L(Largeur element), x(espace entre element)   
% Vecteur X
% [H,L,x]

X = [2/1000 1.5/1000 1/1000 0 0 0];             %point de départ
nb_var = length(X);
dh = 0.001;     %pas pour les dérivés

r = 1;
lambda = [0 0 0];


while r<800000
      
    %MÉTHODE DU GRADIENT
    S = directionS(X, lambda, r, nb_var);
    ii=1;
    while any(abs(S)>0.001)
        %MÉTHODE DE NEWTON
        h = 1;
        dAdh = (LA(X(ii,:)+(h+dh)*S, lambda,r)-LA(X(ii,:)+(h-dh)*S, lambda,r))/(2*dh);
        while abs(dAdh) >0.01
            d2Adh2 = (LA(X(ii,:)+(h+dh)*S, lambda,r)-2*LA(X(ii,:)+(h)*S, lambda,r)+ LA(X(ii,:)+(h-dh)*S, lambda,r))/(dh^2);
            h = h-dAdh/d2Adh2;
            dAdh = (LA(X(ii,:)+(h+dh)*S, lambda,r)-LA(X(ii,:)+(h-dh)*S, lambda,r))/(2*dh);
        end
        X(ii+1,:) = X(ii,:)+h*S;
        ii = ii+1;
        S = directionS(X(ii,:), lambda, r, nb_var);
    end
    lambda = lambda + 2*r.*(contraintes(X));
    r = 2*r;
end
X(end,:)*1000
Fonctionnelle(X)*1000000
function [ S ] = directionS( X, lambda, r, nb_var )
GEN1693_Projet_var
%Calcul de la direction pour la méthode du gradient
    dh = 0.001;
    vdh = zeros(1,nb_var);

    for ii=1:nb_var
        vdh(ii)=dh;
        S(ii)= (LA(X+vdh,lambda,r) - LA(X-vdh,lambda,r))/2*dh; % CALCUL DE LA DIRECTION  
        vdh = zeros(1,nb_var);
    end
%     S(1)= (LA([X(end,1)+dh X(end,2) X(end,3) X(end,4)],lambda,r) - LA([X(end,1)-dh X(end,2) X(end,3) X(end,4)],lambda,r))/2*dh; % CALCUL DE LA DIRECTION  
%     S(2)= (LA([X(end,1) X(end,2)+dh X(end,3) X(end,4)],lambda,r) - LA([X(end,1) X(end,2)-dh X(end,3) X(end,4)],lambda,r))/2*dh; % CALCUL DE LA DIRECTION  
%     S(3)= (LA([X(end,1) X(end,2) X(end,3)+dh X(end,4)],lambda,r) - LA([X(end,1) X(end,2) X(end,3)-dh X(end,4)],lambda,r))/2*dh; % CALCUL DE LA DIRECTION  
%     S(4)= (LA([X(end,1) X(end,2) X(end,3) X(end,4)+dh],lambda,r) - LA([X(end,1) X(end,2) X(end,3) X(end,4)-dh],lambda,r))/2*dh; % CALCUL DE LA DIRECTION  

%     S(1)= (LA([X(1)+dh X(2) X(3) X(4)],lambda,r) - LA([X(1)-dh X(2) X(3) X(4)],lambda,r))/2*dh; % CALCUL DE LA DIRECTION  
%     S(2)= (LA([X(1) X(2)+dh X(3) X(4)],lambda,r) - LA([X(1) X(2)-dh X(3) X(4)],lambda,r))/2*dh; % CALCUL DE LA DIRECTION  
%     S(3)= (LA([X(1) X(2) X(3)+dh X(4)],lambda,r) - LA([X(1) X(2) X(3)-dh X(4)],lambda,r))/2*dh; % CALCUL DE LA DIRECTION  
%     S(4)= (LA([X(1) X(2) X(3) X(4)+dh],lambda,r) - LA([X(1) X(2) X(3) X(4)-dh],lambda,r))/2*dh; % CALCUL DE LA DIRECTION  
end


function [ LAx ] = LA( X, lambda, r )
GEN1693_Projet_var
LAx = Fonctionnelle(X);
g = contraintes(X);

%Assemblage du Lagrange augmenté;
for ii = 1:length(g)
    LAx = LAx + lambda(ii)*g(ii)+r*g(ii)^2;
end
end

function [ F ] = Fonctionnelle( X )
%Fonctionnelle
GEN1693_Projet_var

H = X(end,1);
L = X(end,2);
x = X(end,3);

N = (Ltot - x)/(L + x);
K0 = 2*N*kTE*L^2*kp/(2*N*kTE*L^2 + kp*H);

F = alp^2 *( K*mcp/(((K+2*K0) + (K0*K/(K+K0)) + mcp)) )^2 / (4*pTE*H^2);
end

function [ g ] = contraintes( X )
%contraintes du problème d'optimisation
GEN1693_Projet_var
g(1) = -X(end,3) + 0.5/1000 + X(end,6)^2;
g(2) = -X(end,2) + 1.2/1000 + X(end,5)^2;
g(3) = X(end,1) - 5/1000 + X(end,4)^2;
g=[g(1) g(2) g(3)];
end
