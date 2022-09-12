function u0 = FormRHS(X)
    l = size(X);
    l=l(1);
    u0 = zeros(l^2,1);
    for i = 1:l
        for j = 1:l
            u0(getInd(i,j,l),1) = X(i,j);
        end 
    end
end

function ind = getInd(i,j, m)
    ind = i + (j-1) * m ;
end
