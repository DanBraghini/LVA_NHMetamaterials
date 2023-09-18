function [kL_sem_PB,kL_sem_SB] = function_SEM_EBRing(wv)
%% Created by Danilo Braghini 03.02.2022
%% Spectral element dynamic stiffness matrix for Love's curved beams
global  rho1 rho2 A1 A2 E1 E2 Iz1 Iz2 eta R s01 s02 ne H_fu
N=length(wv);
%% defining grid indexis and matrices sizes
% number of elements on each cell
% ne;
% number of nodes
nn=ne+1;
% number of states on each node
ns=3;% ns=size(kb2,1);     
% number of internal dofs on the internal nodes of the cell
ni=ns*(nn-2); 
% number of degrees of freedom on both nodes of a segment
% nv=2*ns;
% size of cell's dynamic stiffness matrix
nKc = ns*nn; 
% index of measured dof for feedback
n1=ns+1:2*ns;
% index of actuated dof for feedback
n2=2*ns+1:3*ns; 
kL_sem_PB=zeros(2*ns,N);kL_sem_SB=kL_sem_PB;
 % This loop runs through the frequency vector
 for n=1:N 
     w=wv(n);
    %     %% select the damping model
    %     % viscous damping
    %     if dmodel == 'v'
    %         B=rho*c^2;
    %         k_local=sqrt((w^2*rho - 1i*w*eta)/B);
    %     % structural damping (histeresis)
    %     elseif dmodel == 's' 
    %         k_local= w/c;
    %         k_local=k_local*(1-1i*eta/2);
    %     else 
    %         disp('invalidy damping model')
    %         return
    %     endr
    E1=E1*(1+1i*eta);    E2=E2*(1+1i*eta);
    %% build dynamic stiffness matrix on generalized displacement by generalized force
    % Segment 1 = Segment 3
    kr1 = function_SEM_elemEBRing(R,A1,Iz1,rho1,E1,s01,w);
    %% Segment 2  
    kr2 = function_SEM_elemEBRing(R,A2,Iz2,rho2,E2,s02,w);
     %% build cell matrix
     % combining L1 and L2 and L1
    Kc=zeros(nKc,nKc);
    cont=1;
    for i=0:ns:nKc-2*ns
        if mod(cont,2) ~= 0
            Kc(i+1:i+2*ns,i+1:i+2*ns)=Kc(i+1:i+2*ns,i+1:i+2*ns)+kr1;
        else
            Kc(i+1:i+2*ns,i+1:i+2*ns)=Kc(i+1:i+2*ns,i+1:i+2*ns)+kr2;
        end
        cont=cont+1;
    end
    %% adding feedback law
    %Applying shear force V as a function of bending displacement w
    if ne==3
        gain=evalfr(H_fu,1i*w);
        Kc(n2,n1)=Kc(n2,n1)-gain;
    end
    %% condensation
    % lines index
    ls_u=1:ns;
    ls_c=ns+1:ns+ni;
    % ns+ni = nKc-ns
    ls_b=nKc-ns+1:nKc;
    % columns index
    cs_l=1:ns;
    cs_c=ls_c;
    cs_r=ls_b;
    % dividing cells dynamic stiffness matrix Kc in nine block of matrices
    %----------------------------------------------------------%  
    %    Kc = [  kul    |     Kuc    |     Kur   
    %           ------------------------------
    %            kcl    |     Kcc    |     Kcr
    %           ------------------------------
    %            kbl    |     Kbc    |     Kbr ] 
    %----------------------------------------------------------%
    Kul = Kc(ls_u,cs_l);% \in R^(ns x  ns)
    Kuc = Kc(ls_u,cs_c);% \in R^(ns x  ni)
    Kur = Kc(ls_u,cs_r);% \in R^(ns x  ns)
    Kcl = Kc(ls_c,cs_l);% \in R^(ni x  ns)
    Kcc = Kc(ls_c,cs_c);% \in R^(ni x  ni)
    Kcr = Kc(ls_c,cs_r);% \in R^(ni x  ns)
    Kbl = Kc(ls_b,cs_l);% \in R^(ns x  ns) 
    Kbc = Kc(ls_b,cs_c);% \in R^(ns x  ni)
    Kbr = Kc(ls_b,cs_r);% \in R^(ns x  ns)
    Z1=Kcc\Kcl;
    Z2=Kcc\Kcr;

    Kr=[Kul-Kuc*Z1              Kur-Kuc*Z2
          Kbl-Kbc*Z1             Kbr-Kbc*Z2];
    %% Transfer matrix
    % dividing Kr in four block matrices K11,K12,K21,K22 \in R^(ns x ns)
    %----------------------------------------------------------%      
    %     Kr = [  k11    |     K12
    %           ------------------
    %             k21    |    K22]    
    %----------------------------------------------------------%  
    nk=size(Kr,1);
    K11=Kr(1:ns,1:ns);
    K12=Kr(1:ns,ns+1:nk);
    K21=Kr(ns+1:nk,1:ns);
    K22=Kr(ns+1:nk,ns+1:nk);
    Z=K12\K11;
    
    Tc=[-K12\K11         K12\eye(ns)
        K22*Z-K21       -K22/K12];

    %% solving eigenvalue problem
    [~,Lambda]=eig(Tc);
    % floquet multipliers
    %kL(:,n)=log(diag(lambda))/(-im);
    Fm=diag(Lambda);
    Fm=sort(Fm);
    mu=1i*log(Fm);
    %% resultant normalized wavenumber kL
    kL_sem_PB(:,n)=real(mu);
    kL_sem_SB(:,n)=imag(mu);
   
 end
 
end
