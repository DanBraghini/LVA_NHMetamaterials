function u=function_idealfilter(t,Ti,Tf,a)
N=length(t);
for i=1:N
   if t(i) >= Ti && t(i) <= Tf
       u(i)=exp(-1i*t(i)*a);
   else
       u(i)=0;
   end
end
