function output = function_FRF_Ac_v2(dmodel, rho, L1, L2, A1, A2, c, eta, gamma_c, H_pv, ideal_filter, wv,ncell,e, boundary,F)
% by Danilo Braghini 
% last update 08.12.2021
% Using the spectral element dynamic stiffness matrix for Active acoustic
% duct(homogenious duct) with local feedback
%---------------------------------------%
% external load in volume velocity
%---------------------------------------%
N=length(wv);
%% number of spectral elements on each cell
ne = 3;
Lc=L2+2*L1;
Les = Lc/ne;
x = 0:Les:ncell*Lc;
ndofo = ncell*ne+1;

%% Assemblyng dynamic stiffness matrix and solving u(\omega)   
Dg=zeros(ndofo,ndofo,N); Dg_passive=Dg;
if boundary==1
   ndof=ndofo-2;
elseif boundary==2
  ndof=ndofo-1;
elseif boundary==0
  ndof=ndofo;
end
 %% define where the external load is applied
% force in the middle
if e == 'm'
    m = floor(ndof/2);
    if mod(ndof,2) == 0
         F = [zeros(ndof-m-1,1);F;zeros(ndof-m,1)];
    else
        F = [zeros(m,1);F;zeros(m,1)];
    end
% force on the left end
elseif e == 'l'
    F = [F; zeros(ndof-1,1)];
% force on the right end    
elseif e=='r'
    F = [zeros(ndof-1,1);F];
else 
    disp('invalid force position')
    return
end
u=zeros(ndof,N);u_passive=u;u_a=u;
kDg=zeros(1,N);
% This loop runs through the frequency vector
 for n=1:N 
    w=wv(n);
    %% select the damping model
    % viscous damping
    if dmodel == 'v'
        B=rho*c^2;
        k_local=sqrt((w^2*rho - 1i*w*eta)/B);
    % structural damping (histeresis)
    elseif dmodel == 's'  || dmodel == 'n'
        k_local= w/c;
        k_local=k_local*(1-1i*eta/2);
    else 
        disp('invalidy damping model')
        return
    end
    %% build dynamic stiffness matrix
    k11=1+exp(-2*1i*k_local*L1);
        k12= -2*exp(-1i*k_local*L1);
        k21=k12;
        k22=k11;
    Kd1=[ k11   k12
               k21   k22 ].*(A1/(rho*c*(1-exp(-2*1i*k_local*L1))));

    k11=1+exp(-2*1i*k_local*L2);
        k12= -2*exp(-1i*k_local*L2);
        k21=k12;
        k22=k11;

    Kd2=[ k11   k12
              k21   k22 ].*(A2/(rho*c*(1-exp(-2*1i*k_local*L2))));
     % auxiliary variables
     C1=Kd1(2,2)+Kd2(1,1);
    %% adding feedback law
    gain=evalfr(H_pv,1i*w);
    if ideal_filter == 1
        fase=0;
        gain=gain*function_idealfilter(w,2*pi*0.5e3,2*pi*2e3,fase);
    end
     %% adding feedback law
    C2=Kd2(2,1)-gamma_c*gain;
    
    %% combining L1 and L2 and L1 to build cell matrix
    K_cell=[Kd1(1,1)   Kd1(1,2)       0                 0
                  Kd1(2,1)    C1          Kd2(1,2)           0
                   0        C2           C1                Kd1(1,2)                    
                    0        0          Kd1(2,1)          Kd1(2,2)];
        
    K_cell_passive=[Kd1(1,1)   Kd1(1,2)         0                0
                    Kd1(2,1)    C1          Kd2(1,2)             0
                     0         Kd2(2,1)     C1                Kd1(1,2)                    
                     0           0          Kd1(2,1)          Kd1(2,2)];
%% assembling global stiffness matrix
    Dgn=Dg(:,:,n);
    Dgn_p=Dg_passive(:,:,n);
    for i=(1:ne:ndofo-ne)
        Dgn(i:i+ne,i:i+ne)=Dgn(i:i+ne,i:i+ne)+K_cell;
    end
    for i=(1:ne:ndofo-ne)
        Dgn_p(i:i+ne,i:i+ne)=Dgn_p(i:i+ne,i:i+ne)+K_cell_passive;
    end
%% applying boundary conditions
    if boundary==1
        Dgn(1,:)=[];Dgn(:,1)=[];Dgn(end,:)=[];Dgn(:,end)=[];
        Dgn_p(1,:)=[];Dgn_p(:,1)=[];Dgn_p(end,:)=[];Dgn_p(:,end)=[];
    elseif boundary==2
        Dgn(:,1)=Dgn(:,1)+Dgn(:,end); Dgn(1,:)=Dgn(1,:)+Dgn(end,:);
        Dgn(:,end)=[];Dgn(end,:)=[];
        Dgn_p(:,1)=Dgn_p(:,1)+Dgn_p(:,end); Dgn_p(1,:)=Dgn_p(1,:)+Dgn_p(end,:);
        Dgn_p(:,end)=[];Dgn_p(end,:)=[]; 
    end
%% solve the linear system
    u(:,n)= Dgn\F;
    % PS: integrate F if the external load is given in volume
    % acceleration
    u_a(:,n)= 1i*w*Dgn\F;
    kDg(n) =cond(Dgn);
    u_passive(:,n) = Dgn_p\F;
 end
%% useful outputs
     output.kDg=kDg;
     output.ps_v=u;
     output.ps_a=u_a;
     output.xs=x;
     output.ps_passive=u_passive;
end
