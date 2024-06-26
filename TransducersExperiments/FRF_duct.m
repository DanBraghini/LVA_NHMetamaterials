% Last update: 27/11/23 Danilo
% This scirpt uses SEM model for active acoustic duct with local feedback inputs
% dmodel : flag with the chosen damping model
%     's': structural;
%     'v': viscousrho;
%     'n': none;
% eta: damping parameter;
% L0, L1, L2, L3: lengths of the segments of the structure L0-L1-L2-L3;
% A0, A1, A2, A3: cross sectional areas of thestructure;
% kL: value of wavenumber times length. kL \in  [-pb*pi, pb*pi] \subset \R;
% c: sound velocity within the dute (incompressible and perfect/ideal gas model for air);
clear all
close all
clc
[w_1,FRF_1]=readuff('HspeakerOpen.uff');
%[~,FRF_2]=readuff('HspeakerClosed.uff');
ele=load('frf_mic_bk_speaker.mat');
x = linspace(0,2560,1025);                 %Frequência em Hz
%r_inst = abs(inst.FRF_mic_bk_Point3.y_values.values);    %Resposta do microfone de eletreto em Pa/V
FRF_2 = abs(ele.FRF_mic_bk_speaker.y_values.values);   %Resposta do microfone de instrumentação em Pa/V

%[~,FRF_2]=readuff('HspeakerClosed.uff');
[~,FRF_3]=readuff('HmicroOpen.uff');
[~,FRF_4]=readuff('HmicroClosed.uff');
ind1=find(abs(w_1-200)<=1);
ind2=find(abs(w_1-2500)<=2);
ind1x=find(abs(x-200)<=1);
ind2x=find(abs(x-2500)<=2);
fv=w_1(ind1:ind2);
wv=fv*2*pi;
FRF_speakerO=FRF_1(ind1:ind2);
FRF_speakerC=FRF_2(ind1x:ind2x);
FRF_micO=FRF_3(ind1:ind2);
FRF_micC=FRF_4(ind1:ind2);
ncell=1;
%boundary==0: open-open
%boundary==1: opren-closed
boundary=0;
%% acoustic metamaterial set up
dmodel='v';
if dmodel == 'v'
    eta=100;
elseif dmodel == 's'
   eta=0.05;
elseif dmodel == 'n'
    eta=1e-8;
end
 %% number of spectral elements
nse = 3;   
% values for air as an ideal gas at 25ºC
%rho=1.1839;
%c=346.13;
p = 101.325*1e3;                        % Pressão atmosférica [kPa]
T = 273.15+25;                          % Temperatura ambiente a 24ºC [K]
Resp = 287.058;                         % Constante específica dos gases para o ar seco [J/(kg*K)]
rho = p/(Resp*T);                     % Densidade do ar [kg/m^3]
c = sqrt(1.4*p/rho);                    % B=p*gamma(T)
r3=42e-3/2;
A3=pi*r3^2;A0=A3;
r = 33.5e-3/2;
A = pi*r^2;
if boundary==0
    L2 =18.5e-2+0.6*r;
    ndof=nse+1;
    mic=nse;
elseif boundary==1
    L2=18.5e-2;
    ndof=nse+2;
    mic=ndof-2;
end
L0=38.5e-3;
L1=1-18.5e-2;
L3=50e-3;
N=length(wv);
B=rho*c^2;
%% Assemblyng dynamic stiffness matrix and solving u(\omega)   
 %% define where the external load is applied
% force on the left end    
if boundary==0
        F = [1; zeros(ndof-2,1)];
        PV_SEM=zeros(N,ndof-1);PA_SEM=PV_SEM;
 elseif boundary==1
       F = [1; zeros(ndof-1,1)];
       PV_SEM=zeros(N,ndof);PA_SEM=PV_SEM;
 end
% This loop runs through the frequency vector
 for n=1:N 
   Dgn=zeros(ndof,ndof);
    w=wv(n);
    %% select the damping model
    % viscous damping
    if dmodel == 'v'
        % Direct Dynamic Stiffness Matrix Method for  volume velocity
        k_local=sqrt((w^2*rho - 1i*w*eta)/B);
        
        k11=1+exp(-2*1i*k_local*L0);
            k12= -2*exp(-1i*k_local*L0);
            k21=k12;
            k22=k11;
        D0=[ k11   k12
                   k21   k22 ].*(A0/(rho*c*(1-exp(-2*1i*k_local*L0))));
               
        k11=1+exp(-2*1i*k_local*L1);
            k12= -2*exp(-1i*k_local*L1);
            k21=k12;
            k22=k11;
        D1=[ k11   k12
                   k21   k22 ].*(A/(rho*c*(1-exp(-2*1i*k_local*L1))));
               
            k11=1+exp(-2*1i*k_local*L2);
            k12= -2*exp(-1i*k_local*L2);
            k21=k12;
            k22=k11;
        D2=[ k11   k12
                   k21   k22 ].*(A/(rho*c*(1-exp(-2*1i*k_local*L2))));
        
        if boundary==0
               D={D0, D1, D2};
        elseif boundary==1
                k11=1+exp(-2*1i*k_local*L3);
                k12= -2*exp(-1i*k_local*L3);
                k21=k12;
                k22=k11;
                D3=[ k11   k12
                       k21   k22 ].*(A3/(rho*c*(1-exp(-2*1i*k_local*L3))));
                 D={D0, D1, D2, D3};
        end
    % structural damping (histeresis)
    elseif dmodel == 's'  || dmodel == 'n'
        k_local= w/c;
        k_local=k_local*(1-1i*eta/2);
        % Transfer Matrix Method for volume velocity
        Am=[0                                         -1j*k_local*sqrt(B*rho)/A0
              -1j*k_local*A0/sqrt(rho*B)          0];
         T=expm(Am*L0); 
         D0=-[T(1,1)/T(1,2)                                -1/T(1,2)
              T(2,1)-T(2,2)*T(1,1)/T(1,2)         T(2,2)/T(1,2)];
        Am=[0                                         -1j*k_local*sqrt(B*rho)/A
              -1j*k_local*A/sqrt(rho*B)          0];
        T=expm(Am*L1);
        D1=-[T(1,1)/T(1,2)                                -1/T(1,2)
              T(2,1)-T(2,2)*T(1,1)/T(1,2)         T(2,2)/T(1,2)];
        T=expm(Am*L2);
        D2=-[T(1,1)/T(1,2)                                -1/T(1,2)
              T(2,1)-T(2,2)*T(1,1)/T(1,2)         T(2,2)/T(1,2)];
      
          if boundary==0
               D={D0, D1, D2};
        elseif boundary==1
              Am=[0                                         -1j*k_local*sqrt(B*rho)/A3
                  -1j*k_local*A3/sqrt(rho*B)          0];
              T=expm(Am*L3);
              D3=-[T(1,1)/T(1,2)                                -1/T(1,2)
                  T(2,1)-T(2,2)*T(1,1)/T(1,2)         T(2,2)/T(1,2)];
            D={D0, D1, D2, D3};
          end
    else 
        disp('invalidy damping model')
        return
    end
       for i=1:ndof-1
            Dgn(i:i+1,i:i+1)=Dgn(i:i+1,i:i+1)+D{i};
       end
%        Dgn=[D0(1,1)     D0(1,2)                       0                                       0                           0 
%                  D0(2,1)      D0(2,2)+D1(1,1)          D1(1,2)                            0                           0
%                        0               D1(2,1)                 D1(2,2)+D2(1,1)     D2(1,2)                           0
%                        0                       0                     D2(2,1)                  D2(2,2)+D3(1,1)       D3(1,2)
%                      0                       0                           0                           D3(2,1)                  D3(2,2)];
    
    %% applying open-end boundary condition;
    if boundary==0
        Dgn(:,end)=[];Dgn(end,:)=[];
    end
    %% solve the linear system
    %PVO_SEM(n,:)= DgnO\Fo;
    PV_SEM(n,:)= Dgn\F;
    % PS: integrate F if the external load is given in volume acceleration
    %PAO_SEM(n,:)= PVO_SEM(n,:)/(1i*w);
    PA_SEM(n,:)= PV_SEM(n,:)/(1i*w);
 end
%% Computing Transductions
T1_A=FRF_speakerC./PA_SEM(:,2);
T1_V=FRF_speakerC./PV_SEM(:,2);
%% plots
figure
plot(fv/1000,20*log10(abs(PV_SEM(:,mic))),'r','MarkerSize',8,'LineWidth',1)
hold on
plot(fv/1000,20*log10(abs(PA_SEM(:,mic))),'k','MarkerSize',8,'LineWidth',1)
if boundary==0
    plot(fv/1000,79.1+20*log10(abs(FRF_speakerO)),'g')
elseif boundary==1
    plot(x(ind1x:ind2x)/1000,79.1+20*log10(abs(FRF_speakerC)),'g')
end
%plot(fv/1000,20*log10(abs(T1_A)),'*k')
%plot(fv/1000,20*log10(abs(T1_V)),'*r')
% ylabel('$ p / G $ [dB]', 'interpreter', 'latex', 'fontsize', 15)
xlabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
%legend('L','R','L filtered','R Filtered','passive')
xlim([0 fv(end)/1000])
box on
set(gcf, 'Color', 'w');
grid on

% figure
% if boundary==0
%     plot(fv/1000,20*log10(abs(FRF_micO)),'k')
% elseif boundary==1
%     plot(fv/1000,20*log10(abs(FRF_micC)),'r')
% end
% xlabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
% set(gca,'TickLabelInterpreter','Latex','fontsize',15);
% %legend('L','R','L filtered','R Filtered','passive')
% xlim([0 fv(end)/1000])
% box on
% set(gcf, 'Color', 'w');
% grid on
