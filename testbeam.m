% Band gaps in beams with two-geometry cell
clear all
%close all

% Section 1
L1=0.2;
h1=0.02;
b1=0.02;
E1=210e9;
rho1=7800;
eta1=0;
%eta1=0.01;
A1=b1*h1;
I1=b1*h1^3/12;

% Section 2
L2=L1;
%h2=0.005;
h2=h1;
b2=b1;
E2=E1;
rho2=rho1;
eta2=eta1;
A2=b2*h2;
I2=b2*h2^3/12;

wv=2:2:5e3;
N=length(wv);
 %% defining grid indexis and matrices sizes
% number of homogenius segments on each cell
nhs=3;
% number of nodes
nn=nhs+1;
% number of states on each homogenius segment 
ns=2;% ns=size(kb2,1);     
% number of internal dofs on the internal nodes of the cell
ni=ns*(nn-2); 
% number of degrees of freedom on both nodes of a segment
% nv=2*ns;
% size of cell's dynamic stiffness matrix
nKc = ns*nn; 
%% Beam spectral element dynamic stiffness matrix
im=sqrt(-1);
kL=zeros(2*ns,N);
for n=1:N
    w=wv(n);
    E1=E1*(1+1i*eta1);    E2=E2*(1+1i*eta2);
    k_local1=((w^2*rho1*A1)/(E1*I1))^0.25;
    kb1=function_SEM_homogenius_EBBeam(E1,I1,L1,k_local1);
    k_local2=((w^2*rho2*A2)/(E2*I2))^0.25;
    kb2=function_SEM_homogenius_EBBeam(E2,I2,L2,k_local2);
    %K1=Elem_Beam_BE(E1,I1,L1,A1,rho1,eta1,w);
    %K2=Elem_Beam_BE(E2,I2,L2,A2,rho2,eta2,w);
    % Cell
    %% build cell matrix
    % combining L1 and L2 and L1
    Kc=zeros(nKc,nKc);
    c=1;
    for i=0:2:nKc-2*ns
        if mod(c,2) ~= 0
            Kc(i+1:i+2*ns,i+1:i+2*ns)=kb1;
        else
            Kc(i+1:i+2*ns,i+1:i+2*ns)=kb2;
        end
        c=c+1;
    end
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
%     K = [  kul    |     Kuc    |     Kur   
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

    [Phi,lambda]=eig(Tc);
    kL(:,n)=log(diag(lambda))/(-im);
end
kLreal=abs(real(kL));
kLimag=-abs(imag(kL));
%%
figure
plot(wv,kLreal,'k+',wv,kLimag,'r+')
xlabel('\omega [rad/s]')
ylabel('kL(\omega) [rad]')
legend('Real','Imag')
% %%
% figure
% plot(wv,kLreal(1,:),'k+',wv,kLimag(1,:),'r+')
% xlabel('\omega [rad/s]')
% ylabel('kL(\omega) [rad]')
% legend('Real','Imag')
% %%
% figure
% plot(wv,kLreal(2,:),'k+',wv,kLimag(2,:),'r+')
% xlabel('\omega [rad/s]')
% ylabel('kL(\omega) [rad]')
% legend('Real','Imag')
% %%
% figure
% plot(wv,kLreal(3,:),'k+',wv,kLimag(3,:),'r+')
% xlabel('\omega [rad/s]')
% ylabel('kL(\omega) [rad]')
% legend('Real','Imag')
% %%
% figure
% plot(wv,kLreal(4,:),'k+',wv,kLimag(4,:),'r+')
% xlabel('\omega [rad/s]')
% ylabel('kL(\omega) [rad]')
% legend('Real','Imag')
