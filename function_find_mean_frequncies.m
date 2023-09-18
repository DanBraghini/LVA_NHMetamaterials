function output = function_find_mean_frequncies(kL_sem_PB,fv)

tol = 0.001;
ind1_ = find(abs(kL_sem_PB) <= tol);
ind2_ = find(abs(kL_sem_PB-pi) <= tol);
while size(ind1_,2) <= 4 
    tol=10*tol;
    ind1_ = find(abs(kL_sem_PB) <= tol);
end
while size(ind2_,2) <=4 
    tol=10*tol;
    ind2_ = find(abs(kL_sem_PB-pi) <= tol);
end

% Next loops run trough the vectors ind
aux = ind1_(1);
 j=2;
 ind0(1)=ind1_(1);
for i =2:length(ind1_)
    if ind1_(i) ~= aux+1
        ind0(j) = aux;
        ind0(j+1)=ind1_(i);
        j = j+2;
    end
    aux =ind1_(i);
end
nPBs = 2;
aux = ind2_(1);
 j=2;
 indpi(1)=ind2_(1);
for i =2:length(ind2_)
    if ind2_(i) >= aux+10
        indpi(j) = aux;
        indpi(j+1)=ind2_(i);
        j = j+2;
        nPBs=nPBs+1;
    end
    aux =ind2_(i);
end
 output.nPBs =nPBs;
if nPBs >= 1
    for i = 1:indpi(1)
        kL_PBs(1,i) = kL_sem_PB(i);
    end
    if nPBs >= 2
        for i = indpi(2):ind0(3)
            kL_PBs(2,i) = kL_sem_PB(i);
        end
    end
        if nPBs>=3
            for i = ind0(4):indpi(3)
                kL_PBs(3,i) = kL_sem_PB(i);
            end
        end
            if nPBs >= 4
                for i = indpi(4):ind0(5)
                    kL_PBs(4,i) = kL_sem_PB(i);
                end
            end
end
    [~,m1] = min(abs(kL_PBs(1,1:indpi(1))-pi/2));

    if nPBs == 1
        m(1) = m1(1);
    elseif nPBs ==2
         m(1) = m1(1);
         m(2) = indpi(2) + m1;
    elseif nPBs ==3
         m(1) = m1(1);
        m(2) = indpi(2) + m1;
        m(3)= ind0(4) + m1;
    elseif nPBs >= 4
        m(1) = m1(1);
        m(2) = indpi(2) + m1;
        m(3)= ind0(4) + m1;
        m(4) = indpi(4) + m1;
    end
    output.m =m;
end