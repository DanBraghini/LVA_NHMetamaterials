function output = function_buildFEM_PZTRod(varargin)
% function output = function_buildFEM_PZTRod(rho_p, L_s, L_a, Lc, A_s, A_a, e_p,alpha_p,E_p_,C,Gamma_c,ne_cell,ncell,varargin)
% builds the FEM model of the system (PZT-Rod metamaterial) with the
% feddback interaction
% inputs: 
%         ne_cell: number of finite elements per unity cell (multiple of 3 to be divided between the 3 segments)
%         ncell: number of cells that make the structure 
%         rho_p, L_s, L_a, Lc, A_s, A_a, e_p,alpha_pp_: passive material
%         properties
%         C: optimal Rayleigh damping matrix for the structure, computed to
%         make SEM and FEM damping models similar
%         Gamma_c: matrix representing the feedback interactions
%         Boundary conditions opitions:
%            boundary == 0 : clamped-clamped (default)
%            boundary == 1 : open-open
%            boundary == 2 : periodic(infinity system)
%         plotpassive options:
%            0: don't give passive system
%outputs
%       n1: index for the last node of segment 1 on the unity cell (sensor
%       node);
%       n2: index for the last node of segment 2 on the unity cell
%       (actuator node);
%       ndof: number of degrees of freedom after applying boundary
%       conditions;
%% parametric inputs
param=struct(varargin{:});
E_p_=param.Young_undamped;
rho_p=param.density_PZT;
L_s=param.sensor_length;
L_a=param.actuator_length;
A_s=param.sensor_cross_area;
A_a=param.actuator_cross_area;
Gamma_c=param.feedback_matrix_actuator;
K_a=param.passive_matrix_actuator;
%B_s=param.coeficient_sensor;
e_p=param.piezoelectric_constant;
alpha_p=param.dielectric_constant;
Lc=param.cell_length;
%wv=param.frequency_vector;
%F=param.impulse_amplitude;
ncell=param.number_cells;
a=param.non_locality;
%eta=param.damping_coef;
ad=param.damping_FEM;
ne_cell=param.number_FEM_elements_cell;

%%
% inputs:
% fc: frequency of the harmonic excitation F
% ne: number of finite elements per unity cell
% np: number of nodes on the mesh per unity cell
% n1: index for the last node of segment 1 on the unity cell
% n2: index for the last node of segment 2 on the unity cell
% ncell: number of cells that make the structure

%FEM
% number of elements for sensor (S) and actuator (A) PZT materials
ne_a = ne_cell/2;ne_s=ne_a;
ne=ne_cell;
np = ne + 1;
% lengths of each element
Le_a=L_a/ne_a; Le_s = L_s/ne_s;
% grid nodes
x1 = L_s; 
%x2=L_c
% indexing grid
n1 = ne_s + 1; 

M=zeros(np,np);K=M;K_aux=M;

Me_s = (rho_p*A_s*Le_s/6).*[ 2 1
                           1 2] ;

Me_a = (rho_p*A_a*Le_a/6).*[ 2 1
                           1 2];

Ke_p = (A_s/Le_s)*(E_p_ + e_p^2/alpha_p).*[1 -1
                                          -1 1 ];
%% assembling the unity cell
for i=1:ne_s
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + Me_s;
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + Ke_p;
end

for i=n1:ne
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + Me_a;
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + Ke_p;
end

% PZT actuator -> K_a contribution
K(n1,n1) = K(n1,n1) + K_a(1,1);
K(n1,np) = K(n1,np) + K_a(1,2);
K(np,n1) = K(np,n1) + K_a(2,1);
K(np,np) = K(np,np) + K_a(2,2);
if param.plotpassive==1
    K_passive = K;
end
% PZT actuator -> Gamma contribution

%  local-control
% K(n1,1) = K(n1,1) + Gamma_c(1,1);
% K(n1,n1) = K(n1,n1) + Gamma_c(1,2);
% K(np,1) = K(np,1) + Gamma_c(2,1);
% K(np,n1) = K(np,n1) + Gamma_c(2,2);  

%  non-local-control
K_aux(n1,1) = K_aux(n1,1) + Gamma_c(1,1);
K_aux(n1,n1) = K_aux(n1,n1) + Gamma_c(1,2);
K_aux(np,1) = K_aux(np,1) + Gamma_c(2,1);
K_aux(np,n1) = K_aux(np,n1) + Gamma_c(2,2);


%% joining cells to make a structure
ndof = ncell*ne+1;
C = ad(1).*M + ad(2).*K;
Ms=zeros(ndof,ndof);Cs=Ms;Ks=Ms;

for i=1:ne:ndof-ne
    Ms(i:i+ne,i:i+ne) = Ms(i:i+ne,i:i+ne) + M;
    Cs(i:i+ne,i:i+ne) = Cs(i:i+ne,i:i+ne) + C;
    Ks(i:i+ne,i:i+ne) = Ks(i:i+ne,i:i+ne) + K;
end
if param.plotpassive==1
    Ks_passive=Ms;
    for i=1:ne:ndof-ne
        Ks_passive(i:i+ne,i:i+ne) = Ks_passive(i:i+ne,i:i+ne) + K_passive;
    end
end
% Coment this for local-control
% for i=1:ne:ndof-(a+1)*ne
%     Ks( a*ne+i:(a+1)*ne+i , i:i+ne ) = Ks( a*ne+i:(a+1)*ne+i , i:i+ne ) + K_aux;
% end
for i=a*ne+1:ne:ndof-ne
     Ks( i:i+ne,i-a*ne:i-(a-1)*ne) = Ks(  i:i+ne,i-a*ne:i-(a-1)*ne) + K_aux;
end
 
%% Spatial vector x 
xA1 = 0:Le_s:x1;
xB1 = x1+Le_a:Le_a:Lc;

xcell = [xA1 xB1];
x = xcell;
for i = 1:ncell-1
    xaux = xcell + i*Lc + Le_s;
    x = [x xaux];
    x(end) = [];
end
%% applying boundary conditions
if param.boundary==1
    Ms(1,:)=[];Ms(:,1)=[];Ms(end,:)=[];Ms(:,end)=[];
    Ks(1,:)=[];Ks(:,1)=[];Ks(end,:)=[];Ks(:,end)=[];
    Cs(1,:)=[];Cs(:,1)=[];Cs(end,:)=[];Cs(:,end)=[];
    ndof=ndof-2;
elseif param.boundary==2
    Ms(:,1)=Ms(:,1)+Ms(:,end);Ms(:,end)=[];
    Ms(1,:)=Ms(1,:)+Ms(end,:);Ms(end,:)=[]; 
    Ks(:,1)=Ks(:,1)+Ks(:,end); Ks(:,end)=[];
    Ks(1,:)=Ks(1,:)+Ks(end,:);Ks(end,:)=[];
    Cs(:,1)=Cs(:,1)+Cs(:,end);Cs(:,end)=[];
    Cs(1,:)=Cs(1,:)+Cs(end,:);Cs(end,:)=[];
    ndof=ndof-1;
end
%% outputs
output.M=M;
output.C=C;
output.K=K;
output.Ms=Ms;
output.Cs=Cs;
output.Ks=Ks;
output.x=x;
output.ndof=ndof;
output.n1=n1;
output.np=np;
if param.plotpassive==1
    output.Ks_passive=Ks_passive;
end

