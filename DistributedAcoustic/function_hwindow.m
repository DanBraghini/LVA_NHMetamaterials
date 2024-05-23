function u=function_hwindow(t,Ti,Tf)
N=length(t);
for i=1:N
   if t(i) >= Ti && t(i) <= Tf
       u(i)=1;
   else
       u(i)=0;
   end
end
