function [kL_sem_PB,kL_sem_SB] = function_SEM_EBBeam(wv)
%% Created by Danilo Braghini 03.02.2022
%% SEM dynamic stiffness matrix for Euler-Bernoulli Beams
global  rho1 rho2 L1 L2 A1 A2 E1 E2 Iz1 Iz2 eta ne gamma_c H_fu ideal_filter
N=length(wv);
  %% defining grid indexis and matrices sizes
    % number of homogenius segments on each cell
    %ne;
    % number of nodes
    nn=ne+1;
    % number of states on each node 
    ns=2;% ns=size(kb2,1);     
    % number of internal dofs on the internal nodes of the cell
    ni=ns*(nn-2); 
    % number of degrees of freedom on both nodes of a segment
    % nv=2*ns;
    % size of cell's dynamic stiffness matrix
    nKc = ns*nn; 
    % index of measured dof for feedback
    n1=ns+1;
    % index of actuated dof for feedback
    n2=ns*2+1; 
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
    %% defining local wavenumbers
    k_local1=((w^2*rho1*A1)/(E1*Iz1))^0.25;
    k_local2=((w^2*rho2*A2)/(E2*Iz2))^0.25;

    %% build dynamic stiffness matrix on generalized displacement by generalized force
    % Segment 1 = Segment 3
    kb1=function_SEM_homogenius_EBBeam(E1,Iz1,L1,k_local1);

    %% Segment 2  
    kb2=function_SEM_homogenius_EBBeam(E2,Iz2,L2,k_local2);
     %% build cell matrix
     % combining L1 and L2 and L1
    Kc=zeros(nKc,nKc);
    cont=1;
    for i=0:ns:nKc-2*ns
        if mod(cont,2) ~= 0
            Kc(i+1:i+2*ns,i+1:i+2*ns)=Kc(i+1:i+2*ns,i+1:i+2*ns)+kb1;
        else
            Kc(i+1:i+2*ns,i+1:i+2*ns)=Kc(i+1:i+2*ns,i+1:i+2*ns)+kb2;
        end
        cont=cont+1;
    end
    %% adding feedback law
    %Applying shear force V as a function of bending displacement w
    gain=evalfr(H_fu,1i*w);
    if ideal_filter == 1
        fase=0;
        gain=gain*function_idealfilter(w,2*pi*0.5e3,2*pi*2e3,fase);
    end
    Kc(n2,n1)=Kc(n2,n1)-gamma_c*gain;
    %% condensation
    ls_u=1:ns;
    ls_c=ns+1:ns+ni;
    % ns+ni = nKc-ns
    ls_b=nKc-ns+1:nKc;
    cs_l=1:ns;
    cs_c=ls_c;
    cs_r=ls_b;
    % dividing K in nine block of matrices
%----------------------------------------------------------%  
%     K = [  Kul    |     Kuc    |     Kur   
%           ------------------------------
%            Kcl    |     Kcc    |     Kcr
%           ------------------------------
%            Kbl    |     Kbc    |     Kbr ] 
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
