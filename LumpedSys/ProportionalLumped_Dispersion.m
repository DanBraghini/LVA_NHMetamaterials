% Controlled 1D lattice dimer
clear all
close all
clc
%%
% stiffness constant
ke=1; 
% damping factor
b=0;%0.01;
% number of masses per cell
l=1;
% non-locality parameter, a=0 for local control f_k = g_p
% (q_{k-a}-q_{k-a11})
a=0;
% choose the masses
m=zeros(1,l);
m(1)=1;
for j=2:l
    m(j)=m(j-1)+1;
end
wn=sqrt(ke/m(1));
% feedback gaon
gp=-1;
muvec=pi*linspace(-l,l,701);
Omega=zeros(2*l,length(muvec)); 
Modes=zeros(2*l,2*l,length(muvec));
Mj=diag(m);
%------------------------------------------------------------------------------------------------------------------%
%------------------Periodic (Infinity) System : Dispersion Relation--------
%----------------------------------------------------------------------------------------------------------------%
for i=1:length(muvec)
        mu=muvec(i);
        % Cell's stiffness matrix
        Kj=toeplitz([2 -1 zeros(1,l-2)]);
        Kj(1,l)=Kj(1,l)-exp(-1i*mu);Kj(l,1)=Kj(l,1)-exp(1i*mu);
        Kj=ke.*Kj;
        % Cell's viscous damping matrix
        Rj=Kj.*(b/ke);
        % Cell's proportional gain matrix
        if l==1 % analytic solution for single band system
            cmu=1-cos(mu);
            emu=1-exp(1i*mu);
            c1=-2*b/m*cmu*1i;
            c2=gp/m*exp(1i*mu*a)*emu-2*ke/m*cmu;
            lambda=roots([1 c1 c2]);
           % f=find(real(lambda) >= 0);
%             if lambda(1)<=eps
%                 omega(:,i)=[0;0];
%             else

                 Omega(:,i)=sort(lambda,'ComparisonMethod','real');
       %     end
        else % general Dynamic stiffness solution
            Gpj=zeros(l);
            n=fix(a/l); r=rem(a,l);
            ind_pos(1).n=l-r+1;
            ind_pos(1).j=-(n+1);
             ind_neg(1).n=ind_pos(1).n-1;
             ind_neg(1).j=ind_pos(1).j;
            k=1;
            while k<=l
                if ind_pos(k).n<0
                    ind_pos(k).n=l;
                    ind_pos(k).j=ind_pos(k).j-1;          
                 elseif  ind_pos(k).n>l
                     ind_pos(k).n=1;
                     ind_pos(k).j=ind_pos(k).j+1;
                end
                if ind_pos(k).j==0
                       Gpj(k,ind_pos(k).n)=Gpj(k,ind_pos(k).n)+1;
                else
                      Gpj(k,ind_pos(k).n)=Gpj(k,ind_pos(k).n)+exp(1i*ind_pos(k).j*mu);
                end
                 if ind_neg(k).n<0
                    ind_neg(k).n=l;
                    ind_neg(k).j=ind_neg(k).j-1;
                elseif  ind_neg(k).n>l
                     ind_neg(k).n=1;
                     ind_neg(k).j=ind_neg(k).j+1;
                 end
                 if ind_neg(k).j==0
                     Gpj(k,ind_neg(k).n)=Gpj(k,ind_neg(k).n)-1;
                else
                     Gpj(k,ind_neg(k).n)=Gpj(k,ind_neg(k).n)-exp(1i*ind_neg(k).j*mu);
                 end
                  k=k+1;
                  ind_pos(k).n=ind_pos(k-1).n+1;
                  ind_pos(k).j=ind_pos(k-1).j;
                  ind_neg(k).n=ind_pos(k-1).n;
                  ind_neg(k).j=ind_pos(k-1).j;
            end
             Gpj=gp.*Gpj;
%             kc1=gp;kc2=gp;
%     Kc=diag([kc1 kc2])*[(-1)^a1*exp(-1i*mu*a1p)    (-1)^(a1+1)*exp(-1i*mu*a1pp)
%                        (-1)^(a2+1)*exp(-1i*mu*a2p) (-1)^a2*exp(-1i*mu*a2pp)];
%     if norm(Kc-Gpj,'fro')>1e-15 
%         Kc
%         Gpj
%     end
     Kt=Kj-Gpj;%Gpj;
     [V,lambda] = polyeig(Kt,1i.*Rj,-Mj);
        % The Quadratic Eigenvalue problem has 2*np eigenvalues, but we shall take 
        % only those with positive real(omega)
      %  f=find(real(lambda) >= 0);
%        w_aux = lambda(f);
        Omega(:,i) = sort(lambda,'ComparisonMethod','real');
        for ii=1:2*l
            if abs(real(Omega(ii,i)))<=10*eps
                Omega(ii,i)=1i*imag(Omega(ii,i));
            end
            if abs(imag(Omega(ii,i)))<=10*eps
                Omega(ii,i)=real(Omega(ii,i));
            end
        end
        
       %  [V,lambda] = eig(Kt,Mj);
       %  omega(:,i) = sqrt(diag(lambda));
       % [omega(:,i),ind]=sort(omega(:,i));

    % Sort complex numbers according to their real part
    % Modes are the free wave mode shapes
%    [omega(:,i),ind] = sort(lambda(:,i), 'ComparisonMethod','real');
 %   Modes(:,:,i) = V(:,ind); 
        end
end
%---Periodic (Infinit) System : Dispersion Relation------%
%%
figure
Omega=Omega./wn;
for i=2*l:-1:1
    plot3(real(Omega(i,:)),imag(Omega(i,:)),muvec/pi)
    xlabel('$\Omega_r$','interpreter', 'latex','FontSize',15')
    ylabel('$\Omega_i$ ','interpreter', 'latex','FontSize',15')
    zlabel('$\mu/\pi$ ','interpreter', 'latex','FontSize',15')
    box on
    grid on
    set(gcf, 'Color', 'w');
    hold on
    zlim([-1 1])
    xlim([0 max(max(real(Omega)))])
end