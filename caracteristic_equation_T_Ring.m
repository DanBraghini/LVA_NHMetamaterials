%by Danilo Braghini 31.01.2022
% based on SEM formulation of  Danilo Bel, Priscilla Brandao Silva and 
% Prof. Dr. José Roberto de França Arruda 29/08/2013
function f=caracteristic_equation_T_Ring(w)
% Function used to compute the dispersion relation indirectly w(k).
% This function writes down the caracteristic equation of the transfer matrix
% using SEM model for active ring (closed curved beam) according to Euler-Bernoulli
% theory (homogenious) with local feedback.
% As a function of frequency w, corresponding to a fixed wavenumber k,i.e, 
% w(k) the equation is transcendental, since the0 matrix has trigonometric
% or hiperbolic functions (exponentials). Thus, this non-linear equation 
% can be solved by calling the function fsolve on the main code.
% Inputs( given as global variables):
% dmodel : string with the model chosen.
%     's': structural;
%     'v': viscousrho;
%     'n': none;
% L1, L2: lengths of the segments of the unit cell L1-L2-L1;
% A1, A2: cross sectional areas of the unit cell A1-A2-A1;
% kL: value of wavenumber times length. kL \in R | [-pb*pi, pb*pi];
% for pb Brillouin Zones;
% c: sound velocity within the dute;
% eta: damping parameter;
% gamma_c: feedback gain. [k_p k_i k_d] for PID;
% H_pv : transfer fucntion in pressure by volume velocity;
% ideal_filter: 1 if ideal filter is choosen. 0 otherwise.
global  rho1 rho2 L1 L2 A1 A2 kL E1 E2 Iz1 Iz2 gamma_c H_fu ideal_filter

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
nv=2*ns;
% size of cell's dynamic stiffness matrix
nKc = ns*nn; 
% index of measured dof for feedback
n1=ns+1;
% index of actuated dof for feedback
n2=ns*2+1; 

%% defining parameters for local wavenumbers
beta1=((w^2*rho1*A1)/(E1*Iz1))^0.25;
beta2=((w^2*rho2*A2)/(E2*Iz2))^0.25;

%% build dynamic stiffness matrix on generalized displacement by generalized force
% Segment 1 = Segment 3
kb1=function_SEM_homogenius_EBBeam(E1,Iz1,L1,beta1);

%% Segment 2  
kb2=function_SEM_homogenius_EBBeam(E2,Iz2,L2,beta2);
 %% build cell matrix
 % combining L1 and L2 and L1
Kc=zeros(nKc,nKc);
c=1;
for i=0:2:nKc-nv
    if mod(c,2) ~= 0
        Kc(i+1:i+nv,i+1:i+nv)=kb1;
    else
        Kc(i+1:i+nv,i+1:i+nv)=kb2;
    end
    c=c+1;
end
%% adding feedback law
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
%% caracteristic equation
f=det(Tc-exp(-1i*kL).*eye(2*ns));

end


