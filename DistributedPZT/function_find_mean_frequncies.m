function [m1,m2,m3,m4] = function_find_mean_frequncies(kL_sem_PB)

    ind1_ = find(abs(kL_sem_PB) <= 0.01);
    ind2_ = find(abs(kL_sem_PB-pi) <= 0.01);
    
    aux = ind1_(1);
    j=1;
    for i =1:length(ind1_)
        if ind1_(i) ~= aux+1
            ind1(j) = ind1_(i);
            j = j+1;
        end
        aux =ind1_(i);
    end

    aux = ind2_(1);
    j=1;
    for i =1:length(ind2_)
        if ind2_(i) ~= aux+1
            ind2(j) = ind2_(i);
            j = j+1;
        end
        aux =ind2_(i);
    end

    for i = 1:ind2(1)
        kL_PBs(1,i) = kL_sem_PB(i);
    end
    for i = ind2(1):ind1(2)
        kL_PBs(2,i) = kL_sem_PB(i);
    end
    for i = ind1(2):ind2(2)
        kL_PBs(3,i) = kL_sem_PB(i);
    end
    for i = ind2(2):ind1(3)
        kL_PBs(4,i) = kL_sem_PB(i);
    end

    [~,m1] = min(abs( abs(kL_PBs(1,1:ind2(1)))-pi/2 ));
    m1 = m1(1);
    m2 = ind2(1) + m1;
    m3= ind1(2) + m1;
    m4 = ind2(2) + m1;
end