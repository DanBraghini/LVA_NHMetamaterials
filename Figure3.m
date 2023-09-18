%% clear 
clear all 
%close all 
clc 
%% Created on : 15/02/2021 by Danilo Braghini
flim=3.5e3;
dmodel='s';
e='m';
%% Crystal Set up
% acoustic metamaterial
if dmodel == 's'
    eta=0.035;
else
   eta= 10;
end
Lc = 50e-2;
r1 = 2e-2;
A1 = pi*r1^2;
r2 = r1;
A2 = A1;
xs = Lc/4;
xi = xs+Lc/2;
rho =1.225;
c=343;
B=rho*c^2;
L1= xs;
L2 = xi- xs;
%% PID Feedback control law and sensor locality paramater 
kp=1e-7;
ki=1.5e-2;
kd=2e-10;
gamma_c=[0*kp 0*ki kd];
a=0;
%% FEM
% ncell:number of cells that make the structure (multiple of 3 to be divided between the 3 segments)
% number of elements per cell
ne=7*3;
%number of cells on the structure
ncell=18;
% boundary == 0 : free-free (default)
% bundary == 1: fixed-fixed
output = function_buildFEM(rho, L1, L2, Lc, A1, A2, c, eta,ne,ncell,flim,xs,xi,'dmodel',dmodel);
M=output.M;C=output.C;K=output.K;
x=output.x;ndof=output.ndof;
n1=output.n1;n2=output.n2;
% Inserting local feedback interactions 
gamma_c_ph=(B*rho).*gamma_c;
M(n2,n1) = M(n2,n1)-gamma_c_ph(3);
K(n2,n1) = K(n2,n1)-gamma_c_ph(2);
C(n2,n1) = C(n2,n1)-gamma_c_ph(1);
%% joining cells to make a structure
Ms=zeros(ndof,ndof);Cs=Ms;Ks=Ms;

for i=1:ne:ndof-ne
    Ms(i:i+ne,i:i+ne) = Ms(i:i+ne,i:i+ne) + M;
    Cs(i:i+ne,i:i+ne) = Cs(i:i+ne,i:i+ne) + C;
    Ks(i:i+ne,i:i+ne) = Ks(i:i+ne,i:i+ne) + K;
end

%% Transient Response: applied force
% Nc is the number of circles from central frequency fc.
% T2 define the envelope frequency f2, which has Nc circles within.
% Nt is the number of entries on time vector t. The heavside function
% defines a window on half a period of the envelope.
% time discretization period dt depends on N, in such manner that
% the total time of the simulation T can be modified if nedeed.
Nt=4096*2;
% T=0.01; % total analysis time
Ap=1;
fc=1618;
wc = 2*pi*fc;
Nc=15; 
Tc=1/fc;
T2=Nc*Tc*2; 
f2=1/T2; w2=2*pi*f2;
dt=T2/2/4096;
%dt=T2/Nt;
t=0:dt:(Nt-1)*dt;
% t=0:dt:1000e-6;
F=Ap*sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t); 
%  df=1/t(end); 
% fm=1/dt;
%  fv=0:df:fm/2;
% % fN=N*df/2;
% % f=0:df:fN-1;
%  Fw=1/Nt*fft(F);
%   figure
% subplot(1,2,1)
% plot(t,F)
% title('external force dynamics')
% xlabel('time [s]')
% ylabel('F (t)')
% subplot(1,2,2)
% plot(fv/1000,abs(Fw(1:Nt/2)))
% title('external force spectrum')
% xlabel('frequency (kHz)')
% ylabel('|F(\omega)|')
% xlim([0 300])
%% Transient Response: Simualation by integration of SS

Ass=[zeros(ndof,ndof) eye(ndof,ndof)
-inv(Ms)*Ks -inv(Ms)*Cs];

if e == 'm'
    m = floor(ndof/2);
    if mod(ndof,2) == 0
         Fm = [zeros(ndof-m-1,1);1;zeros(ndof-m,1)];
    else
        Fm = [zeros(m,1);1;zeros(m,1)];
    end
% excitation on the left end
elseif e == 'l'
    Fm = [1; zeros(ndof-1,1)];
% excitation on the right end    
elseif e=='r'
    Fm = [zeros(ndof-1,1);1];
else 
    disp('invalid force position')
    return
end

Bss=[zeros(ndof,1)
     (B*rho).*inv(Ms)*Fm];

C1=[eye(ndof,ndof) zeros(ndof,ndof)];
D=0;
sys_x = ss(Ass,Bss,C1,D);

Y=lsim(sys_x,F,t,zeros(2*ndof,1));
%% 
Nplots=100;
DI=round(length(t)/Nplots);
tend=length(t);
figure
for i=1:DI:tend
    plot3(x,(t(i).*10^3)*ones(1,length(x)),Y(i,:),'b');
     hold on
end
xlabel('Length [m]'), ylabel('Time [ms]'), zlabel('Pressure p(x,t)[Pa]')

hold off
set(gcf, 'Color', 'w');
%%
figure
for i=1:1:length(t)
    normm = max(abs(Y(i,:)));
    Y(i,:)=abs(Y(i,:))/normm;
end
mesh(t(100:end).*1e3,x,Y(100:end,:).');
set(gcf, 'Color', 'w');
colormap jet
h=colorbar ;
ylabel(h,'$||p(x,t_i)||_{\infty}$','interpreter', 'latex','FontSize',15)
hold on
ylabel('$x [m]$','interpreter', 'latex','FontSize',15); 
xlabel('$t [ms]$','interpreter', 'latex','FontSize',15');
zlabel('$|p| [Pa]$','interpreter', 'latex','FontSize',15);
view([0,90])
set(gcf, 'Color', 'w');
xlim([0 t(end)*1e3])
box on