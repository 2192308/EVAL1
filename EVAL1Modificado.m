% Author: Mansour Torabi
% Email: smtoraabi@ymail.com
%%  
clc, clear

syms th1 th2  Dth1 Dth2        % permite declarar más de una variable simbólica
syms l1 l2 m1 m2 J1 J2  g t b1 b2   % Declaracion de variables con factores de disipacion

m1v = m1*exp(-t); % Relaciona las masas a usar con una funcion exponencial decreciente respecto al tiempo
m2v = m2*exp(5*t);


%% Kinetic and Potential Energy (Energia potencial y cinetica)
T1 = 1/2*J1*Dth1^2 + 1/2*m1v*(l1/2*Dth1)^2; 
%  J1 es el momento de inercia en el pendulo1;
%  Dth1 es la derivada de theta1
%  m1 es la masa del pendulo1;
%  l1 es la longitud del pendulo1.

Vc2_x = l1*Dth1*cos(th1) + l2/2*(Dth2)*cos(th2);
%  Velocidad en X del pendulo2 al derivar la posicionX del centro de masa.
%  Tiene un componente que depende de la longitud del pendulo1 y del coseno
%  del angulo que este forme.

Vc2_y = l1*Dth1*sin(th1) + l2/2*(Dth2)*sin(th2);
%  Velocidad en Y del pendulo2 al derivar la pisicionY del centro de masa.
%  Tiene un componente que depende de la longitud del pendulo1 y del seno
%  del angulo que este forme.

Vc2 = sqrt(Vc2_x^2 + Vc2_y^2);
%  Se normaliza la velocidad del pendulo2.

T2 = 1/2*J2*(Dth2)^2 + 1/2*m2v*Vc2^2;
%  J2 es el momento de inercia en el pendulo2;
%  Dth2 es la derivada de theta2
%  m2 es la masa del pendulo2;
%  Vc2 al cuadrado es la velocidad del pendulo2 en ambos ejes normalizados. 

T = T1 + T2;
%  Energia cinetica total del sistema.
%  T = 1/2*m1*(Dx1.^2+Dy1.^2)+1/2*m2*(Dx2.^2+Dy2.^2)+1/2*J1*Dth1.^2+1/2*J2*Dth2.^2

V1 = m1v*g*l1/2 * (1-cos(th1));
%  Energia potencial del pendulo1 considerando la gravedad que afecta a Y.
%  m1 es la masa del pendulo1.
%  g es la gravedad de 9.8 m/s2.
%  l1 es la longitud del pendulo1.

V2 = m2v*g*(l1*(1-cos(th1)) + l2/2*(1-cos(th2)));
%  Energia potencial del pendulo2 considerando la gravedad que afecta a Y.
%  m2 es la masa del pendulo1.
%  g es la gravedad de 9.8 m/s2.
%  l1 es la longitud del pendulo1 mas la mitad de la longitud del pendulo2.

V = V1 + V2;
%  Energia potencial de todo el sistema 

D = (1/2)*b1*(Dth1)^2 + (1/2)*b2*(Dth2)^2;% Disipacion

L = T - V; % Lagrangiano 
%%
q  = [th1, th2]; % q hace referencia a los grados de libertad y en que se referencia
Dq = [Dth1, Dth2];% crea un vector simbolico de las velocidades 
tt = linspace(0,5,500); % crea un vector de 0 a 5 con 500 datos igualmente espaciados 
Eq = LagrangeDynamicEqDeriver(L, q, Dq, D); % Llama la funcion de derivacion lagrangiano
[SS, xx] = DynamicEqSolver(Eq, q, Dq, [l1 l2 m1 m2 J1 J2 g b1 b2],...
          [0.5, 0.5, 5, 5, 0.2, 0.2, 9.81, 5, 0], tt, [120,30,0,0]/180*pi); % llama a la funcion que resuelve la ecuacion

% DynamicEqSolver recibe el resultado arrojado por LagrangeDynamicEqDeriver
% recibe a q y a Dq, recibe la lista de variables simbolicas usadas, recibe
% los valores que el ejercicio tenga para cada masa, longitudes etc. 
% recibe el tiempo de estudio tt y las condiciones iniciales propias dadas.
%%             
figure; 
plot(tt,xx(:,1)/pi*180,'r', 'linewidth',2); hold on; plot(tt,xx(:,2)/pi*180,...
    '--b','linewidth',2);

S1 = sprintf('$ \\theta_1$'); 
S2 = sprintf('$ \\theta_2$');
H = legend(S1, S2); 
set(H,'interpreter','latex','fontsize',18,'location','SouthWest');

hx = xlabel('Tiempo (s)');   set(hx, 'fontsize', 18);
hy = ylabel('Angulo (grados)'); set(hy, 'fontsize', 18);
set(gca, 'fontsize', 18);
saveas(gcf, 'Pic/Ex1.png')
title('m1v=m1*exp(-t); m2v=m2*exp(5*t)')
subtitle('b1=5; b2=0[Pa*s]')
%%
Animator1(xx(:,1:2), tt)% llama al animador que ejemplifica matlab como evento
