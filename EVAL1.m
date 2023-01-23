% Author: Mansour Torabi 
% Email: smtoraabi@ymail.com
% Editado por Cristian Sierra, Diego Espinel y Miller Gamba 
clc, clear

syms th1 th2  Dth1 Dth2 D       % permite declarar más de una variable simbólica
syms l1 l2 m1 m2 J1 J2  g t b1 b2   % permitiendo la realizacion de operaciones de manera simbolica y no numerica
m1v = m1*exp(-t); %funcion para la perdida de masa   
m2v = m2*exp(-t);
%% Kinetic and Potential Energy 
T1 = 1/2*J1*Dth1^2 + 1/2*m1v*(l1/2*Dth1)^2; 


Vc2_x = l1*Dth1*cos(th1) + l2/2*(Dth2)*cos(th2);


Vc2_y = l1*Dth1*sin(th1) + l2/2*(Dth2)*sin(th2);


Vc2 = sqrt(Vc2_x^2 + Vc2_y^2);


T2 = 1/2*J2*(Dth2)^2 + 1/2*m2v*Vc2^2;


T = T1 + T2;


V1 = m1v*g*l1/2 * (1-cos(th1));


V2 = m2v*g*(l1*(1-cos(th1)) + l2/2*(1-cos(th2)));


V = V1 + V2;
%  Energia potencial de todo el sistema

D = (1/2)*b1*(Dth1)^2 + (1/2)*b2*(Dth2)^2;% Disipacion

L = T - V; % Lagrangiano 
%%
q  = [th1, th2]; 
Dq = [Dth1, Dth2];
tt = linspace(0,5,500); 
Eq = LagrangeDynamicEqDeriver(L, q, Dq, D); 
[SS, xx] = DynamicEqSolver(Eq, q, Dq, [l1 l2 m1 m2 J1 J2 g b1 b2],...
          [0.5, 0.5, 1, 5, 0.2, 0.5, 9.81, 0.1, 0.2], tt, [120,30,0,0]/180*pi);

%DynamicEqSolver recibe el resultado arrojado por LagrangeDynamicEqDeriver
%recibe a q y a Dq, recibe la lista de variables simbolicas usadas, recibe
%los valores que el ejercicio tenga para cada masa, longitudes etc. 
%recibe el tiempo de estudio tt y las condiciones iniciales propias dadas.           
figure; 
plot(tt,xx(:,1)/pi*180,'r', 'linewidth',2); hold on; plot(tt,xx(:,2)/pi*180,...
    '--b','linewidth',2);

S1 = sprintf('$ \\theta_1$'); 
S2 = sprintf('$ \\theta_2$');
H = legend(S1, S2); 
set(H,'interpreter','latex','fontsize',18,'location','SouthWest');

hx = xlabel('Time (s)');   set(hx, 'fontsize', 18);
hy = ylabel('Angles (grados)'); set(hy, 'fontsize', 18);
set(gca, 'fontsize', 18);
saveas(gcf, 'Pic/Ex1.png')
%%
Animator1(xx(:,1:2), tt)% llama al animador que ejemplifica matlab como evento