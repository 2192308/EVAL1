function Eq = LagrangeDynamicEqDeriver(L, q, Dq, D)
% Author: Mansour Torabi
% Email: smtoraabi@ymail.com

% Editado por Cristian Sierra, Diego Espinel y Miller Gamba  

syms t % variable simbolica t

N = length(q);% N del largo de q (grados de libertad)


%% Calculation of L_q = r.L/r.q and L_Dq = r.L/r.Dq
% Si quiero agregar D entonces seria L_DDq = diff(L, Dq(ii));

L_q = sym(zeros(N,1)); 
L_Dq = sym(zeros(N,1)); 

for ii = 1:N 
   L_q(ii) = diff(L, q(ii));
   L_Dq(ii) = diff(L, Dq(ii)); 
end

%% Calculation of  L_Dq_dt = qd/dt( r_Dq ) 
L_Dq_dt = sym(zeros(N,1)); 

for ii = 1:N
  
    for jj = 1:N
        q_dst = [char(q(jj)), '(t)'];
        Dq_dst = ['diff(', q_dst,',t)'];
        L_Dq(ii)  = subs(L_Dq(ii), {q(jj), Dq(jj)}, {str2sym(q_dst), str2sym(Dq_dst)});  
        
    
    L_Dq_fcn     = symfun(L_Dq(ii), t);
    L_Dq_dt(ii)  = diff(L_Dq_fcn, t);
    
    for jj = 1:N
        
        q_orig = [char(q(jj)), '(t)']; 
        Dq_orig = ['diff(', q_orig,',t)'];
        DDq_orig = ['diff(', q_orig,',t,t)'];
        
        DDq_dst = ['DD',char(q(jj))];
        
  L_Dq_dt(ii)= subs(L_Dq_dt(ii),{str2sym(q_orig), str2sym(Dq_orig), str2sym(DDq_orig)}, ...
                {q(jj), Dq(jj), str2sym(DDq_dst)});
        
    end
end
%% Disipación 
D_Dq = sym(zeros(N,1));

for ii = 1:N
    D_Dq(ii) = diff(D, Dq(ii));
end
%% Lagrange's equations (Second kind) 
Eq = sym(zeros(N,1));

for ii = 1:N
   Eq(ii) = simplify(L_Dq_dt(ii) - L_q(ii) + D_Dq(ii)) ;
end
