function x = function_findx_mesh_NHmetastructure(m, Lc, ne1,ne2,ne3, L1,L2,L3)
    ne=ne1+ne2+ne3;
    Le1=L1/ne1;
    Le2=L2/ne2;
    Le3=L3/ne3; 
    
    cells=floor(m/ne);
    r = mod(m,ne);
    q=r*ne;
    if q <= ne1
        x = cells*Lc+q*Le1;
    else
        if q <=ne1+ne2
            x = cells*Lc+L1+(1-ne1)*Le2;
        else
            x= cells*Lc+L1+L2+(cells-ne1-ne2)*Le3;
        end
    end

end