function output=function_filter(k1,k2,w1,w2,type)
%% low pass filter design
s=tf('s');
if type == 'lp'
    A=10^(-k1/10)-1;
    B=10^(-k2/10)-1;
    n=log10(A/B)/(2*log10(w1/w2));
    n=ceil(n);
    m=1/2/n;
    wc=w2/(B^m);
    i=1;
    k=0;
    while pi*(1+2*k)*m < 0.99*(pi*m+2*pi)
        pk=1i*wc*exp(1i*pi*(1+2*k)*m);
        if real(pk) < 0
             polos(i)=pk;
             i=i+1;
        end
        k=k+1;
    end
     %% designed low pass filter
     Hlp=zpk([],polos,1);
     %% setting static gain to 0 dB
     k=abs(evalfr(Hlp,0));
     Hlp=Hlp*1/k;
     Hf=minreal(Hlp);
     [num,dem]=tfdata(Hf,'v');
     %% eliminating numerical imaginary residues
     num=real(num);
     dem=real(dem);
     %%
     nn=size(num,2);
     nd=size(dem,2);
     Ns=0;Ds=0;
     for i=1:nn 
        Ns=Ns+num(i)*s^(nn-i);
     end
     for i=1:nd
        Ds=Ds+dem(i)*s^(nd-i);
     end
     H=zpk(Ns/Ds);
     Hf=minreal(H);

elseif type == 'bp'
     %% pass band design
     wl=1.1*w1;
     wu=0.9*w2;
     C1=abs((w1^2-wl*wu)/(wl*(wu-wl)));
     C2=((w2^2-wl*wu)/(w2*(wu-wl)));
     w1=1;w2=min(C1,C2);
     A=10^(-k1/10)-1;
     B=10^(-k2/10)-1;
     n=log10(A/B)/(2*log10(w1/w2));
     n=ceil(n);
     m=1/2/n;
     wc=w2/(B^m);
     i=1;
     k=0;
     while pi*(1+2*k)*m < 0.99*(pi*m+2*pi)
        pk=1i*wc*exp(1i*pi*(1+2*k)*m);
        if real(pk) < 0
             polos(i)=pk;
             i=i+1;
        end
        k=k+1;
     end

     H=zpk([],polos,1);
     k=abs(evalfr(H,0));
     H=H*1/k;
     s1=(s^2+wl*wu)/(s*(wu-wl));
     [num,dem]=tfdata(H,'v');
     %% eliminating numerical imaginary residues
     num=real(num);
     dem=real(dem);
     %%
     nn=size(num,2);
     nd=size(dem,2);
     Ns=0;Ds=0;
     for i=1:nn 
        Ns=Ns+num(i)*s1^(nn-i);
     end
     for i=1:nd
        Ds=Ds+dem(i)*s1^(nd-i);
     end
     Hbp=zpk(Ns/Ds);
     Hf=minreal(Hbp);
 elseif type == 'hp'
     %% high pass design
     A=10^(-k1/10)-1;
     B=10^(-k2/10)-1;
     n=log10(A/B)/(2*log10(w1/w2));
     n=ceil(n);
     m=1/2/n;
     wc=w2/(B^m);
     i=1;
     k=0;
     while pi*(1+2*k)*m < 0.99*(pi*m+2*pi)
        pk=1i*wc*exp(1i*pi*(1+2*k)*m);
        if real(pk) < 0
             polos(i)=pk;
             i=i+1;
        end
        k=k+1;
     end
     H=zpk([],polos,1);
     k=abs(evalfr(H,0));
     H=H*1/k;
     s1=w1*w2/s;
     [num,dem]=tfdata(H,'v');
     %% eliminating numerical imaginary residues
     num=real(num);
     dem=real(dem);
     %%
     nn=size(num,2);
     nd=size(dem,2);
     Ns=0;Ds=0;
     for i=1:nn 
        Ns=Ns+num(i)*s1^(nn-i);
     end
     for i=1:nd
        Ds=Ds+dem(i)*s1^(nd-i);
     end
     Hbp=zpk(Ns/Ds);
     Hf=minreal(Hbp);
end
 %% outputs
 output.Hf=Hf;
 output.flporder=n;
end
