function output=function_buildSS(Ms,Cs,Ks,ndof,ncell,n1,n2,e,a,feedback)
    %% inputs

    %% number of elements per cell
    ne_cell=(n2-1)+(n1-1);
    %% Dynamic matrix
    Ass=[zeros(ndof,ndof) eye(ndof,ndof)
        -inv(Ms)*Ks -inv(Ms)*Cs];
    %% External load to state matrix
    % excitation point
    if e == 'm'
        m = floor(ndof/2);
        if mod(ndof,2) == 0
             F = [zeros(ndof-m-1,1);1;zeros(ndof-m,1)];
        else
            F = [zeros(m,1);1;zeros(m,1)];
        end
    % excitation on the left end
    elseif e == 'l'
        F = [1; zeros(ndof-1,1)];
    % excitation on the right end    
    elseif e=='r'
        F = [zeros(ndof-1,1);1];
    else 
        disp('invalid force position')
        return
    end
    B1ss = [zeros(ndof,1)
             inv(Ms)*F];
    %% feedback input to state matrix
    T=zeros(ndof,ncell-a);
    for i=1:size(T,1)
        for j=1:size(T,2)
            if i==n2+(j+a-1)*ne_cell
                T(i,j)=1;
            end
        end
    end

    B2ss = [zeros(ndof,ncell-a)
             inv(Ms)*T];
    %% state to output matrix     
    C1ss=[eye(ndof,ndof) zeros(ndof,ndof)];
    %% state to measured states matrix
    Y=zeros(ncell-a,ndof);
    for i=1:size(Y,1)
        for j=1:size(Y,2)
            if j==n1+(i-1)*ne_cell
                Y(i,j)=1;
            end
        end
    end
    % feedback==1
        C2ss=[Y zeros(ncell-a,ndof)];
  %  elseif feedback==3
 %       C2ss=[zeros(ncell-a,ndof) Y];
  %  end
        
    %% External load to output matrix
    D11ss=zeros(ndof,1);
    %% Feedback input to output matrix
    D12ss=zeros(ndof,ncell-a);
    %% External load to measured states matrix
    D21ss=zeros(ncell-a,1);
    %% outputs
    output.Ass=Ass;
    output.B1ss=B1ss;
    output.B2ss=B2ss;
    output.C1ss=C1ss;
    output.C2ss=C2ss;
    output.D11ss=D11ss;
    output.D12ss=D12ss;
    output.D21ss=D21ss;
end