 function Hf=function_buildfilter(filter, plotfilter,filterinfo)
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
end