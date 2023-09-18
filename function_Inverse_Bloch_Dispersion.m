function  w_TM = function_Inverse_Bloch_Dispersion(fun,options, kLv, wv,pb, a)

global gamma_c
% numerical tolerance 
tol=1e-2;
Nk=length(kLv); 
Nw=length(wv);
w_TM=zeros(pb,Nk);
%% first initial condition 
b=1;
i=1;
m=wv(1);
% saving gamma_c value without spatial delay of non-local feedback
gamma_c_aux =gamma_c;
for n =1:Nk
    kL=kLv(n);
    % for non-local feedback, introduce spatial delay of a cells
    gamma_c = gamma_c_aux.*exp(1i*a*kL);
     w_TM(b,n) = fsolve(fun,wv(m),options);
end
b=b+1;
i=i+1;
%% loop for initial conditions, testing all frequencies on wv
aux=w_TM(b,:);
while i<=Nw && b<=pb
    m=wv(i);
    % loop for wavenumber, computing w(k) for all wavenumbers on kLv
    for n =1:Nk
        kL=kLv(n);
        % for non-local feedback, add spatial delay on gain
        gamma_c = gamma_c_aux.*exp(1i*a*kL);
        % solve w(k) with initial condition wv(m)
        w_TM(b,n) = fsolve(fun,wv(m),options);
    end
    % if the current solution differs on norm, with tolerance tol
    % from the preceding one, a new band was found 
    if norm(aux-w_TM(b,:),2) > tol
        b=b+1;
        aux=w_TM(b,:);
    end   
    i=i+1;
end
gamma_c =gamma_c_aux;
end
