function output=function_initialconditions_setup(pb,fv)
df=fv(2)-fv(1);
if pb>=1
    f1=find(abs(fv-150) <= df/2);
    m=zeros(1,pb);m(1)=f1;
    if pb>=2
        f2=find(abs(fv-500) <= df/2);m(2)=f2;
        if pb>=3
            f3=find(abs(fv-880.3) <= df/2);m(3)=f3;
            if pb>=4
                f4=find(abs(fv-1184) <= df/2);m(4)=f4;
                if pb>=5
                    f5=find(abs(fv-1580) <= df/2);m(5)=f5;
                      if pb>=6
                        f6=find(abs(fv-1860) <= df/2);m(6)=f6;
                        if pb>=7
                            f7=find(abs(fv-2190) <= df/2);m(7)=f7;
                            if pb>=8
                                f8=find(abs(fv-2620) <= df/2);m(8)=f8;
                                if pb>=9
                                     f9=find(abs(fv-2950) <= df/2);m(9)=f9;
                                    if pb>=10
                                        f10=find(abs(fv-3200) <= df/2);m(10)=f10;
                                    end
                                end
                            end
                        end
                      end
                end
            end
        end
    end
end
output=m;