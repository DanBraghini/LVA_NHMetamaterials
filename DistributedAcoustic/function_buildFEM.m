function output = function_buildFEM(rho, L1, L2, Lc, A1, A2, c, eta,ne_cell,ncell,flim,xs,xi,varargin)
% function output = function_buildFEM(rho, L1, L2, Lc, A1, A2, c, eta,ne_cell,ncell,flim,xs,xi,varargin)
% builds the FEM model of the system (acoustic metamaterial) without the
% feddback interaction
% inputs: 
%         ne_cell: number of finite elements per unity cell (multiple of 3 to be divided between the 3 segments)
%         ncell: number of cells that make the structure 
% 
%         Boundary conditions opitions:
%            boundary == 0 : clamped-clamped (default)
%            boundary == 1 : open-open
%            boundary == 2 : periodic(infinity system)
%         Damping model options:
%             viscous: dmodel='v'
%             structural: dmodel='s'
%outputs
%       n1: index for the last node of segment 1 on the unity cell (sensor
%       node);
%       n2: index for the last node of segment 2 on the unity cell
%       (actuator node);
%       ndof: number of degrees of freedom after applying boundary
%       conditions;
%% optional inputs

if nargin > 13
    if nargin == 14
        options = varargin{1};
    else
        options = struct(varargin{:});
    end   
else
    options = [];
end

if ~isfield(options,'boundary')
    options.boundary = 0;
end

if ~isfield(options,'dmodel')
    options.dmodel = 's';
end
dmodel=options.dmodel;
%% FEM
% number of elements
ne_1 = ne_cell/3; ne_2 = ne_1;
%ne==ne_cell
ne = ne_2 + 2*ne_1;
% np: number of nodes on the mesh per unity cell
np = ne + 1;
% lengths of each element
Le_2=L2/ne_2; Le_1 = L1/ne_1;
% indexing grid
n1 = ne_1 +1; n2 = ne_1 + ne_2 +1;
M=zeros(np,np);K=M;
B=rho*c^2;

% mass matrixes
Me1 = rho*A1*Le_1/6*[ 2 1
                      1 2] ;
Me2 = rho*A2*Le_2/6*[ 2 1
                     1 2];  
% viscous damping matrixes   
if dmodel=='v'
    C=M;
    Ce1=eta*A1*Le_1/6*[2 1
                      1 2];                     
    Ce2=eta*A2*Le_2/6*[2 1
                       1 2];  
end
% stiffness matrixes
Ke1 =A1*B/Le_1*[1 -1
                -1 1 ];
Ke2 =A2*B/Le_2*[1 -1
                -1 1 ];
                                                             
                    
%% assembling the unity cell
% this is a loop through elements, not nodes. (each i is an element in a segment, adding a 2x2 matrix) 
for i=1:ne_1
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me1;
    if dmodel=='v'
       C(i:i+1,i:i+1) = C(i:i+1,i:i+1)+ Ce1;
    end
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+ Ke1;
end
for i=n1:ne_1+ne_2
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me2;
    if dmodel=='v'
        C(i:i+1,i:i+1) = C(i:i+1,i:i+1)+ Ce2;
    end
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+ Ke2;
end
for i=n2:ne
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me1;
    if dmodel=='v'
        C(i:i+1,i:i+1) = C(i:i+1,i:i+1)+ Ce1;
    end
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+ Ke1;
end
% structural/histerisis damping model or no damping (eta=0)
if dmodel ~= 'v'
    fm=flim/2;
    C=eta.*K/fm/2/pi;
end

%% joining cells to make a structure
ndof = ncell*ne+1;
Ms=zeros(ndof,ndof);Cs=Ms;Ks=Ms;

for i=1:ne:ndof-ne
    Ms(i:i+ne,i:i+ne) = Ms(i:i+ne,i:i+ne) + M;
    Cs(i:i+ne,i:i+ne) = Cs(i:i+ne,i:i+ne) + C;
    Ks(i:i+ne,i:i+ne) = Ks(i:i+ne,i:i+ne) + K;
end

%% Spatial vector x 
xA1=0:Le_1:xs;
xB1=xs+Le_2:Le_2:xi;
xA2=xi+Le_1:Le_1:Lc;

xcell = [xA1 xB1 xA2];
x = xcell;
for i = 1:ncell-1
    xaux = xcell + i*Lc + Le_1;
    x = [x xaux];
    x(end) = [];
end
%% applying boundary conditions
if options.boundary==1
    Ms(1,:)=[];Ms(:,1)=[];Ms(end,:)=[];Ms(:,end)=[];
    Ks(1,:)=[];Ks(:,1)=[];Ks(end,:)=[];Ks(:,end)=[];
    Cs(1,:)=[];Cs(:,1)=[];Cs(end,:)=[];Cs(:,end)=[];
    ndof=ndof-2;
elseif options.boundary==2
%     Mx=Ms(1:end-1,1:end-1);
%     Mx(1,1)=Mx(1,1)+Ms(end,end); 
%     Mx(1,end)=Mx(1,end)+Ms(end-1,end);
%     Mx(end,1)=Mx(end,1)+Ms(end,end-1);
%     Kx       =Ks(1:end-1,1:end-1);
%     Kx(1,1)=Kx(1,1)+Ks(end,end); 
%     Kx(1,end)=Kx(1,end)+Ks(end-1,end);
%     Kx(end,1)=Kx(end,1)+Ks(end,end-1);
%     Cx       =Cs(1:end-1,1:end-1);
%     Cx(1,1)=Cx(1,1)+Cs(end,end); 
%     Cx(1,end)=Cx(1,end)+Cs(end-1,end);
%     Cx(end,1)=Cx(end,1)+Cs(end,end-1);
%     ndof=ndof-1;
%   Ms = zeros(ndof,ndof); Ks = zeros(ndof,ndof); Cs = zeros(ndof,ndof);
%    Ms = Mx; Ks = Kx; Cs = Cx;
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
output.n2=n2;
output.np=np;

