GEN1693_Projet_var

%%
H=0.5/1000;
L=(1:0.025:2)/1000;
x=(1:0.025:2)/1000;
x=x';

N = (floor((Ltot - x)./(L + x))).^2;
K0 = 2.*N.*kTE.*L.^2*Kp./(2.*N*kTE.*L.^2 + Kp.*H);
deltaT = mcp.*Theta./(((K0.*K)/(K+K0))+mcp);
dT = (K./(K+2.*K0)).* deltaT;

P = N*alp^2*dT.^2/(4.*pTE*(H./L.^2));
V = L.^2*H*N;
F=P./V;

% F = (((floor((Ltot - x)./(L + x))).^2)*alp^2.*((K./(K+2*(2*((floor((Ltot - x)./(L + x))).^2)*kTE.*L.^2*Kp./(2*((floor((Ltot - x)./(L + x))).^2)*kTE.*L.^2 + Kp*H)))).* (mcp*Theta./((((2*((floor((Ltot - x)./(L + x))).^2)*kTE.*L.^2*Kp./(2*((floor((Ltot - x)./(L + x))).^2)*kTE.*L.^2 + Kp*H))*K)/(K+(2*((floor((Ltot - x)./(L + x))).^2)*kTE.*L.^2*Kp./(2*((floor((Ltot - x)./(L + x))).^2)*kTE.*L.^2 + Kp*H))))+mcp))).^2./(4*pTE*(H./L.^2)))./(L.^2*H*(floor((Ltot - x)./(L + x))).^2);

% F=P;
surf(L,x,F,'edgecolor','black');
xlabel('Longueur des modules L (m)')
ylabel('Espace entre modules x (m)')
zlabel('Puissance/volume (W/m)')
title('Puissance diviser par le vilume selon la longueur L et l''espace x')

%%
K = (0:0.1:5);
K0 = (0:0.1:5);
K0 = K0';


deltaT = mcp.*Theta./(((K0.*K)/(K+K0))+mcp);
dT = (K./(K+2.*K0)).* deltaT;
P = N*alp^2*dT.^2/(4.*pTE*(H./L.^2));

surf(K,K0,P)
xlabel('K')
ylabel('K0')
zlabel('Puissance')
title('Impact de K et K0 sur la puissance')






%% 
% 
% H1=1.2/1000;
% L1=1.2/1000;
% x1=(0.5:0.1:2.5)/1000;
% 
% N1 = ((Ltot - x1)./(L1 + x1));
% K01 = 2.*N1*kTE.*L1.^2*Kp./(2.*N1*kTE.*L1.^2 + Kp.*H1);
% 
% F1 = alp^2 *( K*mcp./(((K+2*K01) + (K01*K/(K+K01)) + mcp)) ).^2 / (4*pTE*H1^2);
% 
% H2=(0.5:0.1:2.5)/1000;
% L2=1.2/1000;
% x2=0.5/1000;
% 
% N2 = ((Ltot - x2)./(L2 + x2));
% K02 = 2.*N2*kTE.*L2.^2*Kp./(2.*N2*kTE.*L2.^2 + Kp.*H2);
% 
% F2 = alp^2 *( K*mcp./(((K+2.*K02) + (K02.*K/(K+K02)) + mcp)) ).^2 ./ (4*pTE.*H2.^2);
% 
% H3=1.2/1000;
% L3=(0.5:0.1:2.5)/1000;
% x3=0.5/1000;
% 
% N3 = ((Ltot - x3)./(L1 + x3));
% K03 = 2.*N3*kTE.*L3.^2*Kp./(2.*N3*kTE.*L1.^2 + Kp.*H3);
% 
% F3 = alp^2 *( K*mcp./(((K+2.*K03) + (K03.*K/(K+K03)) + mcp)) ).^2 / (4*pTE.*H3.^2);
% 
% figure(2)
% subplot(1,3,1)
% plot(x1,F1)
% title('F(x)')
% 
% 
% subplot(1,3,2)
% plot(H2,F2)
% title('F(H)')
% 
% subplot(1,3,3)
% plot(L3,F3)
% title('F(L)')




