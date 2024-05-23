function [aM,aK] = function_Calibration_Rayleigh_Damping_non_local(varargin)
% based on the paper: Problems encountered from the use (or misuse) of Rayleigh damping
% To chance where the force is applied, open the function

%% parametric inputs
param=struct(varargin{:});
wv=param.frequency_vector;
eta=param.damping_coef;
[u_SEM,~,~,~,~,~] = function_FRF_PZT(varargin{:});    
u_SEM =u_SEM.';
% figure
% plot(wv/2/pi/1000,20*log10(abs(u_SEM(:,end))),'.-k','MarkerSize',8,'LineWidth',1);hold on;
% plot(wv/2/pi/1000,20*log10(abs(u_SEM(:,1))),'.--r','MarkerSize',8,'LineWidth',1.2)
% xlabel('Frequency [kHz]')
% ylabel('FRF [dB]')
% legend('Left','Right')
% xlim([min(wv/pi/2/1000) max(wv/pi/2/1000)])

[~,LOCS] = findpeaks(abs(u_SEM(:,end)));

% figure
% plot(wv/2/pi/1000,20*log10(abs(u_SEM(:,end))),'.-k','MarkerSize',8,'LineWidth',1);hold on;
% plot(wv(LOCS)/2/pi/1000,20*log10(abs(u_SEM(LOCS,end))),'or','MarkerSize',8,'LineWidth',1.2)
% xlabel('Frequency [kHz]')
% ylabel('FRF [dB]')
% legend('Left','Right')
% xlim([min(wv/pi/2/1000) max(wv/pi/2/1000)])
% title('Passive rod')


wi = wv(LOCS(2));
wj = wv(LOCS(end-5));

aux = (2*eta/(wi+wj))*[wi*wj;1];
alpha = aux(1);
beta = aux(2);

wi = wv(LOCS(2));
for i = 3:length(LOCS)
    R = i-1;
    delta(i) = eta*((1+R-2*sqrt(R))/(1+R+2*sqrt(R)));
    alpha(i) = 2*eta*wi*(2*R/(1+R+2*sqrt(R)));
    beta(i) = 2*eta*(1/wi)*(2/(1+R+2*sqrt(R)));
end
% 
% figure
% plot([1:1:length(LOCS)],alpha,'ob','MarkerSize',5,'LineWidth',1);hold on;grid on;
% xlabel('Number of mode')
% ylabel('alpha')
% 
% figure
% plot([1:1:length(LOCS)],beta,'ob','MarkerSize',5,'LineWidth',1);hold on;grid on;
% xlabel('Number of mode')
% ylabel('beta')
% 
% figure
% plot([1:1:length(LOCS)],100*(delta/eta),'ob','MarkerSize',5,'LineWidth',1);hold on;grid on;
% xlabel('Number of mode')
% ylabel('delta')

aM = alpha(end);
aK = beta(end);

end

