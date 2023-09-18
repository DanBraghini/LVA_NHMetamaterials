clear all
% physical parameters
ke=10; m=1;
gp=-34;
gd=0;
b=1;
% number of degrees of freedom (assuming free-free boundary conditions to build the
% matrices on the standard way)
ndof=9;
M=m.*eye(ndof);
D2_=toeplitz([2 -1 zeros(1,ndof-2)]);
D2_(1,1)=1; D2_(ndof,ndof)=1;
K=ke.*D2_;
C=b.*D2_;
D1_=(diag(-1*ones(1,ndof))+diag(ones(1,ndof-1),-1));
D1_(1,:)=zeros(1,ndof);
Gp=gp.*D1_;
Gd=gd.*D1_;
% imposing fixed-fixed boundary conditions
 ndof_n=ndof-2;
 M(1,:)=[];M(:,1)=[];M(end,:)=[];M(:,end)=[];
 K(1,:)=[];K(:,1)=[];K(end,:)=[];K(:,end)=[];
 C(1,:)=[];C(:,1)=[];C(end,:)=[];C(:,end)=[];
 Gp(1,:)=[];Gp(:,1)=[];Gp(end,:)=[];Gp(:,end)=[]; 
 Gp(1,:)=zeros(1,ndof_n);
 Gd(1,:)=[];Gd(:,1)=[];Gd(end,:)=[];Gd(:,end)=[]; 
 Gd(1,:)=zeros(1,ndof_n);
% closed loop state-space dynamic matrix with x=(q,\dot{q})
 A=[zeros(ndof_n) eye(ndof_n)
      -M\(K+Gp)    -M\(C+Gd)];
%%
% syms b k m g x
% D2_(1,:)=[];D2_(:,1)=[];D2_(end,:)=[];D2_(:,end)=[];
% D1_(1,:)=[];D1_(:,1)=[];D1_(end,:)=[];D1_(:,end)=[];
% X=x*eye(ndof_n)+(b/m+k/m/x).*D2_+g/m/x*D1_;
so=eig(A);
auxm=zeros(size(A));
k=0;
for i=1:length(so)
    auxv=so(i);
    for j=1:length(so)
        if so(j)==auxv && i ~=j
            auxm(i,j)=1;
            k=k+1;
        end
    end
end
[ind_i,ind_j]=find(auxm==1);
% ocorre cruzamento entra os ramos de so(ind_i(k)) e so(ind_j(k)), k \in
% (0,K) com o ganho gp

