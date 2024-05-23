function output = function_buildFEM_xl(m1,m2,k, eta ,ndof,flim, varargin)
% output = function_buildFEM_xl(m1,m2,k, eta ,ndof,flim, varargin)
% builds the FEM model of the system (lumped elastic metamaterial) without the
% feddback interaction
% inputs: 
%         N: number of masses 
% 
%         Boundary conditions opitions:
%            boundary == 0 : free-free (default)
%            boundary == 1 : fixed-fixed
%            boundary == 2 : periodic(infinity system)
%            boundary == 3 : semi-periodic (semi-infinity system)
%outputs

%% optional inputs

if nargin > 6
    if nargin == 7
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


%% FEM
% number of elements
if ndof>1
    mvec=zeros(ndof,1);
    mvec(1:2:ndof,1)=m1; mvec(2:2:ndof,1)=m2;
    Ms=diag(mvec);
    Ks=toeplitz([2*k -k zeros(1,ndof-2)]);
    Ks(1,1)=k;
    Ks(ndof,ndof)=Ks(1,1);
else
    Ms=m1;
    Ks=k;    
end
% structural/histerisis damping model
fm=flim/2;
Cs=eta.*Ks/fm/2/pi;
%% applying boundary conditions
%fixed-fixed
if options.boundary==1
    Ms(1,:)=[];Ms(:,1)=[];Ms(end,:)=[];Ms(:,end)=[];
    Ks(1,:)=[];Ks(:,1)=[];Ks(end,:)=[];Ks(:,end)=[];
    Cs(1,:)=[];Cs(:,1)=[];Cs(end,:)=[];Cs(:,end)=[];
    ndof=ndof-2;
% PBC
elseif options.boundary==2
    Ms(:,1)=Ms(:,1)+Ms(:,end); Ms(1,:)=Ms(1,:)+Ms(end,:);
    Ms(:,end)=[];Ms(end,:)=[]; 
    Ks(:,1)=Ks(:,1)+Ks(:,end);  Ks(1,:)=Ks(1,:)+Ks(end,:);
    Ks(:,end)=[];Ks(end,:)=[];
    Cs(:,1)=Cs(:,1)+Cs(:,end);  Cs(1,:)=Cs(1,:)+Cs(end,:);
    Cs(:,end)=[];Cs(end,:)=[];
    ndof=ndof-1;
end
%% Spatial vector x 
x=1:1:ndof;
%% outputs
output.Ms=Ms;
output.Cs=Cs;
output.Ks=Ks;
output.x=x;
output.ndof=ndof;
end

