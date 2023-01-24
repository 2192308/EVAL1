function Eq = LagrangeDynamicEqDeriver(L, q, Dq, D)
% Author: Mansour Torabi
% Editado por: Cristian Sierra, Miller Gamba y Diego Espinel
% Email: smtoraabi@ymail.com

%%

syms t % variable simbolica t

N = length(q);% N del largo de q (grados de libertad)


%% Calculation of L_q = r.L/r.q and L_Dq = r.L/r.Dq
% Si quiero agregar D entonces seria L_DDq = diff(L, Dq(ii));

L_q = sym(zeros(N,1)); % vector simbolico de ceros de N a 1
L_Dq = sym(zeros(N,1)); % vector simbolico de ceros de N a 1

for ii = 1:N % recorre de 1 a N
   L_q(ii) = diff(L, q(ii));%primera derivada de la funcion L respecto al valor de los angulos: q(ii) osea en la posicion que lleve el for
   L_Dq(ii) = diff(L, Dq(ii)); % primera derivada de la funcion L respecto al valor de la velocidad: Dq(ii) osea en la posicion que lleve el for
end

%% Calculation of  L_Dq_dt = qd/dt( r_Dq ) 
L_Dq_dt = sym(zeros(N,1)); % vector de ceros de N a 1

for ii = 1:N
  
    for jj = 1:N
        q_dst = [char(q(jj)), '(t)'];%crea una variable tipo char con termino (t) en cada posicion de los angulos
        Dq_dst = ['diff(', q_dst,',t)'];%crea una variable tipo char con termino diff(q_dst) "derivada de la velocidad"
        L_Dq(ii)  = subs(L_Dq(ii), {q(jj), Dq(jj)}, {str2sym(q_dst), str2sym(Dq_dst)}); % subs escoge en una expresion y cambia variables por resultados o por otra variable y str2sym evalua lo que tenga dentro del parentesis, evaluando variables simbolicas
        % me dice que de la expresion L_Dq en cada posicion me escoja a la
        % concatenacion de los angulos y las velocidades en cada posicion y
        % me cambie esos valores por la concatenacion de los resultados de
        % cada q_dst y Dq_dst que es la derivada de q y Dq respecto al tiempo
    end
    
    L_Dq_fcn     = symfun(L_Dq(ii), t);% forma de declarar una funcion simbolica del resultado de L_Dq respecto al tiempo
    L_Dq_dt(ii)  = diff(L_Dq_fcn, t);% deriva la funcion simbolica respecto t y la guarda en cada posicion
    
    for jj = 1:N
        
        q_orig = [char(q(jj)), '(t)']; %crea una variable tipo char con termino (t) en cada posicion de los angulos
        Dq_orig = ['diff(', q_orig,',t)'];%crea una variable tipo char con termino diff(q_dst) "derivada de la velocidad" respecto a t
        DDq_orig = ['diff(', q_orig,',t,t)'];
        
        DDq_dst = ['DD',char(q(jj))];
        
        L_Dq_dt(ii)   = subs(L_Dq_dt(ii), {str2sym(q_orig), str2sym(Dq_orig), str2sym(DDq_orig)}, ...
                        {q(jj), Dq(jj), str2sym(DDq_dst)});
        
    end
end
%% Disipacion
D_Dq = sym(zeros(N,1));

for ii = 1:N
    D_Dq(ii) = diff(D, Dq(ii));
end
%% Lagrange's equations (Second kind) 
Eq = sym(zeros(N,1));

for ii = 1:N
   Eq(ii) = simplify(L_Dq_dt(ii) - L_q(ii) + D_Dq(ii)) ;% simplifica la funcion L_Dq_dt para obtener una funcion simplificada sin alterar los terminos
end
