           % Gpsym=zeros(l);
           syms mu
           a=4;l=3;
           Gpsym=zeros(l).*mu;
            n=fix(a/l); r=rem(a,l);
            ind_pos(1).n=l-r+1;
            ind_pos(1).j=-(n+1);
             ind_neg(1).n=ind_pos(1).n-1;
             ind_neg(1).j=ind_pos(1).j;
            k=1;
            while k<=l
                if ind_pos(k).n<0
                    ind_pos(k).n=l;
                    ind_pos(k).j=ind_pos(k).j-1;          
                 elseif  ind_pos(k).n>l
                     ind_pos(k).n=1;
                     ind_pos(k).j=ind_pos(k).j+1;
                end
                if ind_pos(k).j==0
                       Gpsym(k,ind_pos(k).n)=Gpsym(k,ind_pos(k).n)+1;
                else
                      Gpsym(k,ind_pos(k).n)=Gpsym(k,ind_pos(k).n)+exp(1i*ind_pos(k).j*mu);
                end
                 if ind_neg(k).n<0
                    ind_neg(k).n=l;
                    ind_neg(k).j=ind_neg(k).j-1;
                elseif  ind_neg(k).n>l
                     ind_neg(k).n=1;
                     ind_neg(k).j=ind_neg(k).j+1;
                 end
                 if ind_neg(k).j==0
                     Gpsym(k,ind_neg(k).n)=Gpsym(k,ind_neg(k).n)-1;
                else
                     Gpsym(k,ind_neg(k).n)=Gpsym(k,ind_neg(k).n)-exp(1i*ind_neg(k).j*mu);
                 end
                  k=k+1;
                  ind_pos(k).n=ind_pos(k-1).n+1;
                  ind_pos(k).j=ind_pos(k-1).j;
                  ind_neg(k).n=ind_pos(k-1).n;
                  ind_neg(k).j=ind_pos(k-1).j;
            end
            Gpsym