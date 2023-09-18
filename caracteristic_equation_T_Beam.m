%by Danilo Braghini 22.02.2022
% based on SEM formulation
function f=caracteristic_equation_T_Beam(w)
% Function used to compute the dispersion relation indirectly w(k).
% This function writes down the caracteristic equation of the transfer matrix
% using SEM model for active bea) according to Euler-Bernoulli or
% Timoshenko theories, with feedback law giben by Hfu.
% As a function of frequency w, corresponding to a fixed wavenumber k,i.e, 
% w(k) the caracteristic equation is transcendental, since the matrix has trigonometric
% or hiperbolic functions (exponentials). Thus, this non-linear equation 
% can be solved by calling the function fsolve on the main code, with 
% a suitable initial guess.
% Inputs( given as global variables):
% L1, L2: lengths of the segments of the unit cell L1-L2-L1 [m];
% A1, A2: cross sectional areas of the unit cell A1-A2-A1 [m^2];
% Iz1, Iz2: moments of inertia of the  unit cell Iz1-Iz2-Iz1 [m^4];
% kL: value of wavenumber times length. kL \in R | [-pb*pi, pb*pi] 
% when pb is the number of Brillouin Zones;
% E1, E2: Young modulus for the segments of the unit cell E1-E2-E1 [Pa];
% eta: viscoelastic damping parameter;
% Kap: x-sectional constant (Kap=5/6 for rectangular x-section)
% G : Shear modulus [Pa]
% gamma_c: feedback gain. [k_p k_i k_d] for PID;
% H_fu : transfer fucntion in general displacement by general force;
% ideal_filter: 1 if ideal filter is choosen. 0 otherwise.
%-------------------------------------------------------------------------------
%                       SEM mesh 
%-------------------------------------------------------------------------------
% each element is a two node spectral element
% Dynamic stiffness matrix of each element is of the form:
%--------------------------------------------------------
% [f_l = Kc [xs_l 
%  f_r]      xs_r]
%--------------------------------------------------------
% where f(x,t) = [V(x,t) M(x,t)]^T contais the generalized forces
% V = shear internal load
% M = bending moment 
% x(x,t) = [w(x,t) \phi(x,t)]^T represents the state vector at each fixed x an t
%, in this case the generalized displacements:
% w = bending displacement
% \phi = bending angle
% The subscripts _l and _r denote each one of the nodes of the element,
% fixing the space position x. Thus xs_l, xs_r \in R^(2 x 1)
global  rho1 rho2 L1 L2 A1 A2 E1 E2 Iz1 Iz2 eta ne gamma_c H_fu ...
    theory Kap1 Kap2 G1 G2 kL
%% defining grid indexis and matrices sizes
% number of homogenius segments/elements on each cell
% ne
% number of nodes
nn=ne+1;
% number of states on each node
ns=2;% ns=size(kb2,1);     
% number of internal dofs on the internal nodes of the cell
ni=ns*(nn-2); 
% size of cell's dynamic stiffness matrix
nKc = ns*nn; 
% index of measured dof for feedback
n1=ns+1;
% index of actuated dof for feedback
n2=ns*2+1; 
E1=E1*(1+1i*eta);    E2=E2*(1+1i*eta);
%% build dynamic stiffness matrix on generalized displacement by generalized force
% Segment 1 = Segment 3
if theory == 'EB'
    k_local1=((w^2*rho1*A1)/(E1*Iz1))^0.25;
    kb1=function_SEM_homogenius_EBBeam(E1,Iz1,L1,k_local1);
elseif theory == 'T'
    kb1=function_SEM_homogenius_TBeam(w,rho1,A1,E1,Iz1,Kap1,L1,G1);
else
    disp('invalid theory')
end
%% Segment 2  
if theory == 'EB'
    k_local2=((w^2*rho2*A2)/(E2*Iz2))^0.25;
    kb2=function_SEM_homogenius_EBBeam(E2,Iz2,L2,k_local2);
elseif theory == 'T'
    kb2=function_SEM_homogenius_TBeam(w,rho2,A2,E2,Iz2,Kap2,L2,G2);
else
    disp('invalid theory')
end

%% build cell matrix
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
% Applying shear force V as a function of bending displacement w
if ne==3
    gain=evalfr(H_fu,1i*w);
    Kc(n2,n1)=Kc(n2,n1)-gamma_c*gain;
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
%     Kc = [ kul    |     Kuc    |     Kur   
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
