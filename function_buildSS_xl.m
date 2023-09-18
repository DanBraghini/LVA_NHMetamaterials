function output=function_buildSS_xl(Ms,Cs,Ks,e)
% output=function_buildSS_xl(Ms,Cs,Ks,e)
%   builds the State space lumped model based on FEM of the system 
%(lumped elastic metamaterial) without the feddback interaction
%% inputs
%   (Ms,Cs,Js) : FEM discretization of the passive system P 
%       e      : string with excitation point choosen. options are¨middle,
%       left or right ends of the onedimensional structure 
%% number of dof's per cell =1
N=size(Ms,1);
nu=N;nz=N;ny=N;
%% Dynamic matrix
Ass=[zeros(N,N) eye(N)
    -inv(Ms)*Ks -inv(Ms)*Cs];
%% External load to state matrix
% excitation in the middle
if e == 'm'
    m = floor(N/2);
    if mod(N,2) == 0
         F = [zeros(N-m-1,1);1;zeros(N-m,1)];
    else
        F = [zeros(m,1);1;zeros(m,1)];
    end
% excitation on the left end
elseif e == 'l'
    F = [1; zeros(N-1,1)];
% excitation on the right end    
elseif e=='r'
    F = [zeros(N-1,1);1];
else 
    disp('invalid force position')
    return
end
B1ss = [zeros(N,1)
         inv(Ms)*F];
%% feedback input to state matrix
T=[zeros(N-nu,nu)
    eye(nu)];

B2ss = [zeros(N,nu)
         inv(Ms)*T];
%% state to output matrix     
C1ss=[eye(N,N) zeros(N,N)];
%% state to measured states matrix
if ny>1
    Y=diag(-1.*ones(1,ny))+diag(ones(1,ny-1),1);
    Y=[Y zeros(ny,1)];
    Y(ny,ny+1)=1;
else
    Y=[-1 1];
end
C2ss=[Y zeros(ny,2*N-size(Y,2))];
%% External load to output matrix
D11ss=zeros(N,1);
%% Feedback input to output matrix
D12ss=zeros(N,nu);
%% External load to measured states matrix
D21ss=zeros(ny,1);
%% outputs
output.Ass=Ass;
output.B1ss=B1ss;
output.B2ss=B2ss;
output.C1ss=C1ss;
output.C2ss=C2ss;
output.D11ss=D11ss;
output.D12ss=D12ss;
output.D21ss=D21ss;
end