%% clear 
clear all 
close all 
clc 
%echo on
%% lumped-parameter system set up
ke=10; m1=1; m2=1;flim=3.5e0;df=2e-4;ndof=9;
%gp=-37:0.01:6.05;
% gd=(-15:0.01:-5).*1e0;
 %gp=-33:1:-10;
 %gp=-gp;
%gd=-1000:1:-900;
gd=-11:0.01:0;
%gp=0*gd;
gp=-10.*ones(size(gd));
b=[0 0.01 0.05 0.08 0.1 0.12 0.15 0.2:0.1:1 1.5 2];
% angular frequency vector limit
wv=2*pi*(1e-4:df:flim);
fv=wv/2/pi;
r_max=zeros(size(b));g_max=r_max;g_min=r_max;area_max=r_max;ind_a=r_max;gp_amax=r_max; 
r_ga=r_max;norm2_max=r_max;ind_g=r_max;gp_h2max=r_max;
for j=1:length(b)
    % number of elements
         if ndof>1
             mvec=zeros(ndof,1);
             mvec(1:2:ndof,1)=m1; mvec(2:2:ndof,1)=m2;
             Ms=diag(mvec);
             D2_=toeplitz([2 -1 zeros(1,ndof-2)]);
             D2_(1,1)=1; D2_(ndof,ndof)=1;
             Ks=ke.*D2_;
            %viscous damping model
             Cs=b(j).*D2_;
         else
              disp('Warning: only one mass isnt enough to have feedback' )
         end
    %     %fixed-fixed
    %    if boundary==1
    ndof_n=ndof-2;
     Ms(1,:)=[];Ms(:,1)=[];Ms(end,:)=[];Ms(:,end)=[];
     Ks(1,:)=[];Ks(:,1)=[];Ks(end,:)=[];Ks(:,end)=[];
     Cs(1,:)=[];Cs(:,1)=[];Cs(end,:)=[];Cs(:,end)=[];
    so=zeros(2*ndof_n,length(gp));
    q1m=zeros(length(gp),length(wv));q1p=q1m;qnm=q1m;qnp=q1m;
    k=0;
    for i=1:length(gp)
        %gamma_c=[gp(i) 0 0];
       % output=function_CL4MSS(ndof,m1,m2,k,b,flim,boundary,e,gamma_c);
        %% Computing FRF for each gain
         % Feedback Matrix
          D1_=(diag(-1*ones(1,ndof))+diag(ones(1,ndof-1),-1));
          D1_(1,:)=zeros(1,ndof);
          Gsp=gp(i).*D1_;
          Gsd=gd(i).*D1_;
         Gsp(1,:)=[];Gsp(:,1)=[];Gsp(end,:)=[];Gsp(:,end)=[]; 
         Gsp(1,:)=zeros(1,ndof_n);
         Gsd(1,:)=[];Gsd(:,1)=[];Gsd(end,:)=[];Gsd(:,end)=[]; 
         Gsd(1,:)=zeros(1,ndof_n);
 
     %% define where the external load is applied
    %  % F=1 for FRF
      F=1;
     % force in the middle
         m = floor(ndof_n/2);
         if mod(ndof_n,2) == 0
              F = [zeros(ndof_n-m-1,1);F;zeros(ndof_n-m,1)];
         else
             F = [zeros(m,1);F;zeros(m,1)];
         end
        %% closed-loop state-space 4-matrix representation SIMO for FRFs
        nw=1;nz=ndof_n;
        % dynamic matrix
        A=[zeros(ndof_n) eye(ndof_n)
               -Ms\(Ks+Gsp)    -Ms\(Cs+Gsd)];
        % disturbance to state matrix
        B=[zeros(ndof_n,1)
                 Ms\F];
        % state to measued output matrix
        C=[eye(ndof_n)  zeros(ndof_n)];
        % disturbance to measured output matrix
        D=zeros(nz,nw);
        sys_simo=ss(A,B,C,D); 
        %sys_simo=output.CLsys;
        %[A,B,C,D]=ssdata(sys_simo);
        P=tf(sys_simo);
        % Taking BODE diagrams of first and last DOFs
        P1=P(1,1);Pn=P(ndof_n,1);
        % BODE diagram of tf for output q(1)
            [mag,phase]=bode(P1,wv);
            q1m(i,:)=vec(mag);
            q1p(i,:)=vec(phase);
            % BODE diagram of tf for output q(n)
            [mag,phase]=bode(Pn,wv);
            qnm(i,:)=vec(mag);
            qnp(i,:)=vec(phase);
      
                %% Computing Eigenspectum for each gain (free-vibration problem)
                % s=-i \omega
                %((-i\omega)^2 M -i\omega C +K+G)q(\omega) = 0
            %     if b ~=0
            %         [V,so(i,:)]=polyeig(Ks+Gs,Cs,Ms);
            %     else
            %         [V,L]=eig(Ks+Gs,-Ms);
            %         so(i,:)=[sqrt(diag(L));-sqrt(diag(L))];
            %     end
             so=eig(A);
             r=max(real(so));
             if r<-1e-4
                    %% state-space 4-matrix representation SISO for H2 norm
                    nw=1;nz=1;
                    C=[1 zeros(1,ndof_n-2) -1 zeros(1,size(A,2)-ndof_n)];
                    D=zeros(nz,nw);
                    sys_siso{i}.statespace=ss(A,B,C,D);
                     g_stable(k+1)=gd(i);
                     r_stable(k+1)=r;
                     ind_r(k+1)=i;
                     norm2(k+1)=norm(sys_siso{i}.statespace,2);
                     dw=wv(2)-wv(1);
                     area(k+1)=2*(trapz(q1m(i,:))-trapz(qnm(i,:)))*df;
                    % first time inside this condition happens for the first
                     % stable g
                     if k==0
                         ind_g0=i;
                     end
                     k=k+1;
                    
              end
    end
    if k>0
        r_max(j)=max(r_stable);
        g_max(j)=max(g_stable);
        g_min(j)=min(g_stable);
        [area_max(j),ind_a(j)]=max(abs(area));
        if gp(2)~=gp(1)
            g=gp;
        else
            g=gd;
        end  
        g_amax(j)=g(ind_g0-1+ind_a(j));
        r_ga(j)=r_stable(ind_g0-1+ind_a(j));
        [norm2_max(j),ind_g(j)]=max(norm2);
        g_h2max(j)=g(ind_g0-1+ind_g(j));
    end
end
save('G:\Outros computadores\My MacBook Air\MATLAB\NHLumped_sys\TablePD10.mat');

%     figure
%     scatter(gp,abs(area),'k','LineWidth',1.5)
%     hold on
%    % scatter(gp,area,'b','LineWidth',1.5)
%     %scatter(gp,sqrt(abs(area_rms)),'c','LineWidth',1.5)
%     scatter(gp,norm2_matlab,'r','LineWidth',1.5);
%     set(gcf, 'Color', 'w');
%     box on
%     grid on
%     xlabel('$g_p$', 'interpreter', 'latex', 'fontsize', 15')
%     %legend('$|area|$','$area$','$H_2$ norm','interpreter', 'latex')
%     legend('$|area|$','$H_2$ norm','interpreter', 'latex')
%     set(gca,'TickLabelInterpreter','Latex','fontsize',15);
