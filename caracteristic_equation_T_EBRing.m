%by Danilo Braghini 31.01.2022
% based on SEM formulation of  Danilo Bel, Priscilla Brandao Silva and 
% Prof. Dr. José Roberto de França Arruda 29/08/2013
function f=caracteristic_equation_T_EBRing(w)
% Function used to compute the dispersion relation indirectly w(k).
% This function writes down the caracteristic equation of the transfer matrix
% using SEM model for active rings (closed curved beams) according to Euler-Bernoulli
% theory with feedback law giben by Hfu.
% As a function of frequency w, corresponding to a fixed wavenumber k,i.e, 
% w(k) the caracteristic equation is transcendental, since the matrix has trigonometric
% or hiperbolic functions (exponentials). Thus, this non-linear equation 
% can be solved by calling the function fsolve on the main code, with 
% a suitable initial guess.
% Inputs( given as global variables):
% rho1 rho2: mass densities of the segments of the unit cell rho1-rho2-rho3
% [kg/m^3];
% R: ring's radius [m];
% A1, A2: cross sectional areas of the unit cell A1-A2-A1 [m^2];
% Iz1, Iz2: moments of inertia of the  unit cell Iz1-Iz2-Iz1 [m^4];
% s01 s02: curved lengths corresponding to the segments [m];
% kL: value of wavenumber times length. kL \in R | [-pb*pi, pb*pi] 
% when pb is the number of Brillouin Zones;
% E1, E2: Young modulus for the segments of the unit cell E1-E2-E1 [Pa];
% eta: viscoelastic damping parameter;
% gamma_c: feedback gain. [k_p k_i k_d] for PID;
% H_fu : transfer fucntion in general displacement by general force;
% ideal_filter: 1 if ideal filter is choosen. 0 otherwise.
% SEM mesh: 
% each element is a two node spectral element
% Dynamic stiffness matrix of each element is of the form:
%--------------------------------------------------------
% [f_l = Kc [xs_l 
%  f_r]      xs_r]
%--------------------------------------------------------
% where f(x,t) = [N(x,t) V(x,t) M(x,t)]^T contais the generalized forces
% N = longitudinal internal load [N];
% V = shear internal load [N];
% M = bending moment [Nm]; 
% x(x,t) = [u(x,t) w(x,t) \phi(x,t)]^T represents the state vector at each fixed x an t
%, in this case the generalized displacements:
% u = displacement on the tangencial direction of the curved beam [m];
% w = bending displacement, on the perpendicular direction of the curved
% beam [m];
% \phi = bending angle [rad];
% Obs: The subscripts _l and _r denote each one of the nodes of the element,
% fixing the space position x. Thus xs_l, xs_r \in R^(2 x 1)
global  rho1 rho2 R A1 A2 kL E1 E2 Iz1 Iz2 s01 s02 H_fu ne
%% defining grid indexis and matrices sizes
% number of homogenius segments on each cell
%ne
% number of nodes
nn=ne+1;
% number of states on each homogenius segment 
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
%% build dynamic stiffness matrix on generalized displacement by generalized force
% Segment 1 = Segment 3
kr1 = function_SEM_homogenius_EBRing(R,A1,Iz1,rho1,E1,s01,w);
%% Segment 2  
kr2 = function_SEM_homogenius_EBRing(R,A2,Iz2,rho2,E2,s02,w);
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


