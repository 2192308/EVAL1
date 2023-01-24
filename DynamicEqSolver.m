function [SS, xx] = DynamicEqSolver(Eq, q, Dq, ParamList, ParamVal, tspan, InitCnd)
% Author: Mansour Torabi
% Editado por: Cristian Sierra, Miller Gamba y Diego Espinel
% Email: smtoraabi@ymail.com

%% [1.1]: Convert Eq To State-Space Form:

N = length(Eq); % N del largo de Eq

DDq = sym(zeros(1, N));
for ii = 1:N
    DDq(ii) = sym(['DD', char(q(ii))]); %variable simbolica DD de caracter tipo char de los angulos
end
%

% AA * X = BB;

AA = jacobian(Eq, DDq); % saca el jacobiano de Eq respecto a DDq teniendo en cuenta que el vector DDq no tiene relevancia como se exprese en filas o columnas
BB = -simplify(Eq - AA*DDq.'); % simplifica la relacion de las funciones y le asigna un negativo

DDQQ   = sym(zeros(N, 1));
DET_AA = det(AA);% calcula el determinante de la matriz AA que es el resultado del jacobiano

for ii = 1:N   
    AAn       = AA; % iguala a la matriz AA
    AAn(:,ii) = BB;% asigna BB que es la simplificacion negativa de la expresion de antes y se la asigna a cada posicion de AAn hasta el largo de N
    DDQQ(ii)  = simplify(det(AAn) / DET_AA); % simplifica la division entre determinantes y se la asigna a cada posicion de DDQQ hasta el largo de N
end

%% [1.2]: State Space formation - Final Step

SS = sym(zeros(N, 1));

for ii = 1:N
   SS (ii) = Dq(ii); % iguala SS a cada velocidad en cada posicion
   SS (ii + N) = DDQQ(ii);  % iguala la SS en cada posicion mas N a cada posicion de DDQQ
end

%% [1.3]: Change variables from q to x

Q = [q, Dq]; % vectoriza q y Dq en uno solo
X = sym('x_',[1 2*N]); % crea un vecto de "x_" de 1 a dos veces el largo de N
SS = subs(SS, Q, X); % cambia de la exprecion SS los terminos Q por el termino X

%% [2.1] Solving ODEs
syms t
% Preparation of SS Eq for ODE Solver: Creating Anonymous Fcn

SS_0 = subs(SS, ParamList,ParamVal); %cambia de SS los terminos ParamList por ParamVal1
SS_ode0 = matlabFunction(SS_0,'vars',{X, t}); % genera una función anónima de MATLAB a partir del objeto sym SS_0. Las variables libres de SS_0 se convierten en las entradas de la función resultante manejar SS_ode0.
% el vars pide a la funcion anonima unos valores de analisis y {X,t] asigna
% unas condiciones iniciales de la ecuacion

SS_ode  = @(t, x)SS_ode0(x(1:2*N)',t);
[ts, xx] = ode45(SS_ode, tspan, InitCnd);
