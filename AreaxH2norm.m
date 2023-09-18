na=1000;
sys_ss=cell(6,na);norm2=sys_ss;P=sys_ss;area=sys_ss;area_2=sys_ss;err=sys_ss;err2=sys_ss;
flim=300;
df=2e-2;
wv=2*pi*(0:df:flim);
for i=1:6
    for j=1:na
        k=1;
      while k
            sys_ss{i,j}=rss(i+1);
            D=sys_ss{i,j}.D;
            norm2{i,j}=norm(sys_ss{i,j},2);
            norm_inf=norm(sys_ss{i,j},inf);
            %test properness and reasonable norms
            if D == 0 && norm2{i,j}> 5 &&  norm2{i,j}< 70 &&  norm_inf> 5 &&  norm_inf< 70
                k=0;
            else
                continue
            end
      end
            P{i,j}=tf(sys_ss{i,j});
            [mag,phase]=bode(P{i,j},wv);
            Q=vec(mag);
            %area_2{i,j}=sqrt(2*trapz(Q.^2)*df);
            area{i,j}=2*trapz(Q)*df;
            err{i,j}=norm2{i,j}*(sqrt(2*2*pi*flim))-area{i,j};
            %err2{i,j}=norm2{i,j}-area_2{i,j};
            % if norm2{i,j}==inf
           %     err{i,j}=0;err2{i,j}=0;
           % end
    end
end
%%
x=cell2mat(err);
ind=find(abs(x)>=1);
x2=cell2mat(err2);
ind2=find(abs(x2)>=1);
%err>0 => norm>area
ind_upp=find(x>0);
%err<0 => norm<area
ind_low=find(x<0);
ind_inf=find(abs(x)==inf);
P_low=size(ind_low,1)/(size(ind_low,1)+size(ind_upp,1)-size(ind_inf,1));
%%
save('G:\Outros computadores\My MacBook Air\MATLAB\NHLumped_sys\AreaxH2Norm.mat');
