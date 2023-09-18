% UNICAMP - FEM - DMC
% Laboratório de Vibroacústica
%
% Danilo Beli 
% Priscilla Brandao Silva
% Prof. Dr. José Roberto de França Arruda
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Flexible Ring
%%%%%%%%%%%%%%%%%%%%%%%%%
%
%clc
clear all
%close all
%
im = sqrt(-1);          % imaginary number
%

% Ring parameters %%%%%%%%%%%%%%
% R        = 0.050;         % radius [m]
% bh      = 0.010;         % width [m]
% hh      = 0.001;         % height [m]
% rhoh   = 1500;         % mass density [kg/m3]
% Eh       = 5*10^9;     % elastic modulus [Pa]
% etah    = 0.001;         % structural damping
R=0.25;
bh=0.15;
hh=0.002;
rhoh=7200;
Eh=220e9;
etah=0.001;
etah=0;
Ehd = Eh*(1+im*etah);
%
% Spectral Element Mesh %%%%%%%%%%%
Nuc = 4;                      % number of structures or unit cells
Ngdl = 3*Nuc;              % number of degrees of freedom
dtheta = 2*pi/Nuc;       % angular length of the unit cell
theta = 0:dtheta:2*pi-dtheta;   % vector theta
s0 = R*dtheta;               % circumferential length of the unit cell
%
% External Force [u1, w1, phi1, u2, w2, phi2, ...] 
% (1 - longitudinal, 2 - bending, 3 - moment)
Fg     = zeros(Ngdl,1);
Fg(1) = 0; 
Fg(2) = 1; 
Fg(3) = 0;
%
% Auxiliar %%%%%%%%%%%%%%%%%%%
flim=1e3;%in kHz
v_omega=2*pi*(2:2:flim);      % frequência de análise da resposta [rad/s] 0.2
%
    for ct2 = 1:1:length(v_omega)
        omega = v_omega(ct2);
        %
        % Elements
        [D_elem_ring,gamma] = D_elem_flexiblering2(R,bh,hh,rhoh,Ehd,s0,omega);
        %
        % Transfer Matrix -------------------------------------------------
        DLL = D_elem_ring(1:3,1:3);
        DLR = D_elem_ring(1:3,4:6);
        DRL = D_elem_ring(4:6,1:3);
        DRR = D_elem_ring(4:6,4:6);
        invDLR = (DLR\eye(3));
        %T = ([-invDLR*DLL -invDLR;DRL-DRR*invDLR*DLL -DRR*invDLR]); % Transfer Matrix
        T = ([-invDLR*DLL invDLR;-DRL+DRR*invDLR*DLL -DRR*invDLR]); % Transfer Matrix
        [PHI,MI] = eig(T); % eigenvalues and eigenvectors
        [MI,axS] = sort(diag(MI)); PHI = PHI(:,axS);    
        vK(:,ct2) = log((MI))/(-im); % wavenumbers
        %
        % Forced Response %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dg = zeros(Ngdl,Ngdl); 
        for ct1 = 1:1:Nuc-1 % assembly
            ind1=3*(ct1-1)+1; ind2=3*(ct1)+3;
            Dg(ind1:ind2,ind1:ind2) = Dg(ind1:ind2,ind1:ind2)+D_elem_ring;
        end
        % end-end: closing the ring
        Dg(Ngdl-2:Ngdl,Ngdl-2:Ngdl) = Dg(Ngdl-2:Ngdl,Ngdl-2:Ngdl)+D_elem_ring(1:3,1:3);
        Dg(1:3,1:3) = Dg(1:3,1:3)+D_elem_ring(4:6,4:6);
        Dg(1:3,Ngdl-2:Ngdl) = Dg(1:3,Ngdl-2:Ngdl)+D_elem_ring(4:6,1:3);
        Dg(Ngdl-2:Ngdl,1:3) = Dg(Ngdl-2:Ngdl,1:3)+D_elem_ring(1:3,4:6);
        %
        % Displacements - Linear System %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        warning off all;
        Ug(:,ct2) = Dg\Fg;
        warning on all;
    end
%
fsize = 11; fname = 'Times New Roman';
height = 8; width = 6;
ax_p = '-'; ax_LineWidth = 1; ax_MarkerSize = 4;
%
cL = sqrt(real(Eh)/rhoh);
%% 
% figure(1),  plot(v_omega/(2*pi)/cL,real(vK),'.'), box on, 
% xlim([min(v_omega/(2*pi)/cL) max(v_omega/(2*pi)/cL)])
% xlabel('$\overline{\omega}$','interpreter','latex','FontSize',fsize,'FontName',fname), 
% ylabel('$\Re(k)$','interpreter','latex','FontSize',fsize,'FontName',fname), 
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperType','A4','PaperPosition', [0.0 0.0 1 1],'Color',[1 1 1]);
% set(gca,'FontSize',fsize,'FontName',fname,'Box','on','linewidth',1,'YScale','linear');
% %%
% figure(2),  plot(v_omega/(2*pi)/cL,imag(vK),'.'), box on, 
% xlim([min(v_omega/(2*pi)/cL) max(v_omega/(2*pi)/cL)])
% xlabel('$\overline{\omega}$','interpreter','latex','FontSize',fsize,'FontName',fname), 
% ylabel('$\Im(k)$','interpreter','latex','FontSize',fsize,'FontName',fname), 
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperType','A4','PaperPosition', [0.0 0.0 1 1],'Color',[1 1 1]);
% set(gca,'FontSize',fsize,'FontName',fname,'Box','on','linewidth',1,'YScale','linear');
% %%
% figure(3),  plot(v_omega/(2*pi)/cL,20*log10(abs(Ug(1,:)))), box on, 
% xlim([min(v_omega/(2*pi)/cL) max(v_omega/(2*pi)/cL)])
% xlabel('$\overline{\omega}$','interpreter','latex','FontSize',fsize,'FontName',fname), 
% ylabel('$u (dB)$','interpreter','latex','FontSize',fsize,'FontName',fname), 
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperType','A4','PaperPosition', [0.0 0.0 1 1],'Color',[1 1 1]);
% set(gca,'FontSize',fsize,'FontName',fname,'Box','on','linewidth',1,'YScale','linear');
% %
% figure(4),  plot(v_omega/(2*pi)/cL,20*log10(abs(Ug(2,:)))), box on, 
% xlim([min(v_omega/(2*pi)/cL) max(v_omega/(2*pi)/cL)])
% xlabel('$\overline{\omega}$','interpreter','latex','FontSize',fsize,'FontName',fname), 
% ylabel('$u (dB)$','interpreter','latex','FontSize',fsize,'FontName',fname), 
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperType','A4','PaperPosition', [0.0 0.0 1 1],'Color',[1 1 1]);
% set(gca,'FontSize',fsize,'FontName',fname,'Box','on','linewidth',1,'YScale','linear');

%%
figure
plot(v_omega/(2*pi*1000),real(vK)/pi,'k+',v_omega/(2*pi*1000),imag(vK)/pi,'r+');
xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$k L_c/ \pi$', 'interpreter', 'latex', 'fontsize', 15)
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
box on
grid on
set(gcf, 'Color', 'w');

