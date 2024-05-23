function Teff =  function_getMicTransd
% This function gives a effective constant transducer model 
% based on experimental data. It uses average value of the measured FRF
% within the frequency range of interest for our purposes. Uncoment the Figure 
% plots and make manual changes if you need to change this interval.
% The script also passes the FFF from the measured V/Pa (eletreto microphone)
% The measures was taken from the duct on LVA laboratory, at FEM-UNICAMP (Campinas/Brasil)
% 16/07/2022
    load dadosmicrofone.mat
    % figure
    % plot(xeletreto,20*log10(abs(yeletreto)),xinst,20*log10(abs(yinst/13.9e-3))+18.5);
    % legend('eletreto','instrumentacao')
    % axis([0 3000 -inf inf])
    % grid on

    % experimentally evaluated eletreto mic sensitivity
    %s=1.652e-3;
    s=1.5723e-3;
    Teff=s;
end