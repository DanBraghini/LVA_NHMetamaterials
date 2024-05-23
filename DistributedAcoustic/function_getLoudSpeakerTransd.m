function Teff = function_getLoudSpeakerTransd
% This function gives a effective constant transducer model 
% based on experimental data. It uses average value of the measured FRF
% within the frequency range of interest for our purposes. Uncoment the Figure 
% plots and make manual changes if you need to change this interval.
% The script also passes the FRF from the measured mm/s/V (laser
% dopler interferometer)to m^3/s^2/V by multiplying by the average cross
% sectional area of the Loudspeaker's moving part. This average was
% measured from the duct on LVA laboratory, at FEM-UNICAMP (Campinas/Brasil)
% 16/07/2022
    load dadosinterferometro.mat
%    figure
%    subplot(211)
  %  FRF1=1i.*(FRF1).*(2*pi*fFRF1);
%     plot(fFRF1,20*log10(abs(FRF1)),'k',fFRF5,20*log10(abs(FRF5)),'r')
%     grid on
%     legend('1 mm/s/V','5 mm/s/V')
%     xlabel('Frequency: Hz')
%     ylabel('FRF: dB ref mm/s/V')
%     
%     subplot(212)
%     plot(fcoh1,coh1,'k',fcoh5,coh5,'r')
%     grid on
%     legend('1 mm/s/V','5 mm/s/V')
%     xlabel('Frequency: Hz')
%     ylabel('Coherence')

    % The FRF experimentally obtained is in velocity per voltage, thus it is
    % the integral of acceleration per voltage.Take the derivative
    % in order to get the FRF in acceleration per voltage.
    FRF1=1i.*(FRF1).*(2*pi*fFRF1);
    % multiply by the effective cross section area to get
    % volume acceleration per voltage
    d=1e-3*(17.8+26.65)/2;
    FRF1 = FRF1.*(1e-3*pi*d^2/4);
    %figure
     %plot(fFRF1,20*log10(abs(FRF1)),'k')
    % find the first frequency for which coherense becomes close to 1
    ind= find(abs(coh1-1) <=1e-2);ind1=ind(1);
    % impose the frequency for which the range of interest ends
    fend=2000;
    % Take the effective transduction as the mean value
    Teff=(abs(FRF1(ind1))+abs(FRF1(fFRF1 == fend) ) )/2;
end