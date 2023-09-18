function function_metasetup(filter, plotfilter,filterinfo)
%%
if nargin < 4
    if filter == 1
        filterinfo.w1=1800*2*pi;
        filterinfo.w2=2000*2*pi;
        %lp= 'low pass', bp = 'band pass', 'hp' ='high pass'
        filterinfo.type='hp';
    elseif filter == 2
         filterinfo.wf=1e3*2*pi;
    end
end
global dmodel rho L1 L2 Lc A1 A2 c eta gamma_c H_pa kp ki kii kd kdd Nd
%% acoustic metamaterial set up
if dmodel == 's'
    eta=0.001;
elseif dmodel == 'v'
   eta= 100;
elseif dmodel == 'n'
    eta=1e-8;
end
rho=1.225;
c=343;
Lc = 50e-2;
r1 = 2e-2;
A1 = pi*r1^2;
r2 = r1;
A2 = pi*r2^2;
xs = Lc/4;
xi = xs+Lc/2;
L1= xs;
L2 = xi- xs;

%% setting analog filters
% no filter
if filter == 0
    Hf=1;
elseif filter == 1% butterworth filter
    out=function_filter(-3,-10,filterinfo.w1,filterinfo.w2,filterinfo.type);
    Hf=out.Hf;
elseif filter == 2% RC filter 
    wf=filterinfo.wf;
    Hf = s/wf/(1+s/wf);
elseif filter == 3
%------ filter built for the experiment---------------
     Hf=Filtro_ver2_completo_v1;
end
if plotfilter==1
    s1=1i*wv;
    FRFf=zeros(1,length(wv));
    for j=1:length(wv)
         [num,dem]=tfdata(Hf,'v');
         num=real(num);
         dem=real(dem);
         nn=size(num,2);
         nd=size(dem,2);
         Ns=0;Ds=0;
         for i=1:nn 
            Ns=Ns+num(i)*s1(j)^(nn-i);
         end
         for i=1:nd
            Ds=Ds+dem(i)*s1(j)^(nd-i);
         end
         FRFf(j)=Ns/Ds;
    end
    figure
    plot(wv/2/pi/1000,20*log10(abs(FRFf)))
     xlabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
     ylabel('$[dB]$', 'interpreter', 'latex', 'fontsize', 15)
     figure
    plot(wv/2/pi/1000,180*angle(FRFf)/pi)
     xlabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
     ylabel('$degrees$', 'interpreter', 'latex', 'fontsize', 15)
    box on
    set(gcf, 'Color', 'w');
    grid on
end
    
    %% setting feedback law
    s=tf('s');
    % P-I-D-DD in pressure by volume acceleration
    H_pa=[ 1; 1/s; Nd/(1+Nd/s); (Nd/(1+Nd/s))^2; (1/s)^2];
    %% feedback and filter  
    H_pa=H_pa.*Hf;
    %% Transducers: 
    % T1 = sensor, T2 = actuator
    T1 = function_getMicTransd;
    T2 =  function_getLoudSpeakerTransd;
    %T1=1;T2=1;
    H_pa=T2.*H_pa;
    H_pa=H_pa.*T1;
    %% electric gains
    gamma_c=[kp ki kd kdd kii]./(T1*T2);
end