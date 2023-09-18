function son=function_SortMAC(so)
    son=so;
    for i=1:size(so,1)
        for j=1:size(so,2)-1
            rig=son(i,j);
            r_k=zeros(1,size(so,1));
            k=1;
            for k=1:size(so,1)           
                     rkgg=son(k,j+1);
                     r_k(k)=rig'*rkgg/norm(rig,2)/norm(rkgg,2);
            end
             [~,ind_r]=max(r_k);
             aux=son(i,j+1);
             son(i,j+1)=son(i,ind_r);
             son(i,ind_r)=aux;
        end
    end
end