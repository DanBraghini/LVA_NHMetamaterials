function kb=function_SEM_homogenius_EBBeam(E,Iz,L,beta)
% builds dynamic stiffness matriz for a euler bernoulli beam homogenius spectral
% element (SEM), for a given frequency (corresponding to k_local(w))
xsi=beta*L;
im=1i;
    k11=(-1+im)*(xsi^3*exp(2*xsi)-xsi^3*exp(-2*im*xsi)+im*xsi^3*exp((2-2*im)*xsi)...
        -im*xsi^3)/(exp(-2*im*xsi)+exp((2-2*im)*xsi)-4*exp((1-im)*xsi)+1+exp(2*xsi));

    k12=im*(-xsi^3+xsi^3*exp(2*xsi)-exp((2-2*im)*xsi)*xsi^3+xsi^3*exp(-2*im*xsi))...
        /(beta*exp(-2*im*xsi)+beta*exp((2-2*im)*xsi)-4*beta*exp((1-im)*xsi)+beta+beta*exp(2*xsi));

    k13=(-2*im*xsi^3*exp(xsi)-2*xsi^3*exp(-im*xsi)+2*xsi^3*exp((2-im)*xsi)+2*im*xsi^3*...
        exp((1-2*im)*xsi))/((exp(-2*im*xsi)+exp((2-2*im)*xsi)-4*exp((1-im)*xsi)+1+exp(2*xsi)));

    k14=2*xsi^3*(exp(xsi)-exp((2-im)*xsi)-exp(-im*xsi)+exp((1-2*im)*xsi))/(beta*(exp(-2*im*xsi)...
        +exp((2-2*im)*xsi)-4*exp((1-im)*xsi)+1+exp(2*xsi)));



    k21=k12;

    k22=(1+im)*(xsi^3*exp(2*xsi)-im*xsi^3*exp((2-2*im)*xsi)-xsi^3*exp(-2*im*xsi)+im*xsi^3)...
        /(beta^2*exp(-2*im*xsi)+beta^2*exp((2-2*im)*xsi)-4*beta^2*exp((1-im)*xsi)+beta^2+beta^2*exp(2*xsi));

    k23=-2*xsi^3*(exp(xsi)-exp((2-im)*xsi)-exp(-im*xsi)+exp((1-2*im)*xsi))/(beta*((exp(-2*im*xsi)...
        +exp((2-2*im)*xsi)-4*exp((1-im)*xsi)+1+exp(2*xsi))));

    k24=(-2*im*xsi^3*exp(xsi)+2*xsi^3*exp(-im*xsi)-2*xsi^3*exp((2-im)*xsi)+2*im*xsi^3*exp((1-2*im)...
        *xsi))/(beta^2*exp(-2*im*xsi)+beta^2*exp((2-2*im)*xsi)-4*beta^2*exp((1-im)*xsi)+beta^2+beta^2*exp(2*xsi));


    k31=k13;

    k32=k23;

    k33=(1-im)*(xsi^3*exp(-2*im*xsi)-xsi^3*exp(2*xsi)-im*xsi^3*exp((2-2*im)*xsi)+im*xsi^3)/...
        (exp(-2*im*xsi)+exp((2-2*im)*xsi)-4*exp((1-im)*xsi)+1+exp(2*xsi));

    k34=-im*(-xsi^3+xsi^3*exp(2*xsi)-exp((2-2*im)*xsi)*xsi^3+xsi^3*exp(-2*im*xsi))/...
        (beta*exp(-2*im*xsi)+beta*exp((2-2*im)*xsi)-4*beta*exp((1-im)*xsi)+beta+beta*exp(2*xsi));


    k41=k14;

    k42=k24;

    k43=k34;

    k44=(1+im)*(xsi^3*exp(2*xsi)-im*xsi^3*exp((2-2*im)*xsi)-xsi^3*exp(-2*im*xsi)+im*xsi^3)/...
        (beta^2*exp(-2*im*xsi)+beta^2*exp((2-2*im)*xsi)-4*beta^2*exp((1-im)*xsi)+beta^2+beta^2*exp(2*xsi));


    kb=E*Iz/L^3*[ k11  k12  k13  k14
                      k21  k22  k23  k24
                      k31  k32  k33  k34
                      k41  k42  k43  k44 ];
end
                  