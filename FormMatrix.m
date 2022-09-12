
function mat = FormMatrix(u, alpha)
    m = sqrt( length(u));
    beta = 10^(-6);
    A = zeros(m^2, m^2);
    h = 1/(m+1);
    
    for i = 1:m
        for j = 1:m
            Uij=u(getInd(i,j,m));
        
            if i == 1 && j == 1
                Uimj=0;
                Uijm=0;
                Uijp= u (getInd(i,j+1, m));
                Uipj = u(getInd(i+1,j,m));
                Uipjm =0;
                Uimjp = 0;

                AW = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uimjp-Uimj)/h)^2+beta) )  );
                AE = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uipj-Uipjm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AS = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipjm-Uijm)/h)^2 +((Uij-Uijm)/h)^2+beta) )  );
                AN = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) ) +  1/(2*sqrt ( ((Uijp-Uimjp)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AC = 1 - AW - AE - AS - AN;

                
                A( getInd(i,j,m) , getInd(i,j,m) ) = AC;
                A(getInd(i,j,m), getInd(i+1,j,m)) = AE;
                A (getInd(i,j,m), getInd(i,j+1,m)) = AN;
                
            elseif i == m && j == 1
                Uimj=u(getInd(i-1,j,m));
                Uijm=0;
                Uijp= u (getInd(i,(j+1),m));
                Uipj = 0;
                Uipjm =0;
                Uimjp = u(getInd(i-1,j+1,m));

                AW = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uimjp-Uimj)/h)^2+beta) )  );
                AE = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uipj-Uipjm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AS = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipjm-Uijm)/h)^2 +((Uij-Uijm)/h)^2+beta) )  );
                AN = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) ) +  1/(2*sqrt ( ((Uijp-Uimjp)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AC = 1 - AW - AE - AS - AN;



                A( getInd(i,j,m) , getInd(i,j,m) ) = AC;
                A( getInd(i,j,m) , getInd(i-1,j,m) ) = AW;
                A (getInd(i,j,m), getInd(i,j+1,m)) = AN;


            elseif i ==1 && j == m 
   
                Uimj=0;
                Uijm= u(getInd(i,j-1,m));
                Uijp= 0;
                Uipj = u(getInd(i+1,j,m));
                Uipjm = u(getInd(i+1,j-1,m));
                Uimjp = 0;
                AW = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uimjp-Uimj)/h)^2+beta) )  );
                AE = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uipj-Uipjm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AS = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipjm-Uijm)/h)^2 +((Uij-Uijm)/h)^2+beta) )  );
                AN = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) ) +  1/(2*sqrt ( ((Uijp-Uimjp)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AC = 1 - AW - AE - AS - AN;



                A( getInd(i,j,m) , getInd(i,j,m) ) = AC;
                A(getInd(i,j,m), getInd(i+1,j,m)) = AE;
                A(getInd(i,j,m), getInd(i,j-1,m)) = AS;
            
            elseif i == m && j == m
                Uimj= u(getInd(i-1,j,m));
                Uijm=u(getInd(i,j-1,m));
                Uijp= 0;
                Uipj = 0;
                Uipjm = 0;
                Uimjp = 0;
                AW = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uimjp-Uimj)/h)^2+beta) )  );
                AE = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uipj-Uipjm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AS = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipjm-Uijm)/h)^2 +((Uij-Uijm)/h)^2+beta) )  );
                AN = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) ) +  1/(2*sqrt ( ((Uijp-Uimjp)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AC = 1 - AW - AE - AS - AN;


                A( getInd(i,j,m) , getInd(i,j,m) ) = AC;
                A( getInd(i,j,m) , getInd(i-1,j,m) ) = AW;
                A(getInd(i,j,m), getInd(i,j-1,m)) = AS;

            elseif i == 1 
                Uimj=0;
                Uijm=u(getInd(i,j-1,m));
                Uijp= u (getInd(i,(j+1),m));
                Uipj = u(getInd(i+1,j,m));
                Uipjm =u(getInd(i+1,j-1,m));
                Uimjp = 0;
                AW = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uimjp-Uimj)/h)^2+beta) )  );
                AE = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uipj-Uipjm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AS = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipjm-Uijm)/h)^2 +((Uij-Uijm)/h)^2+beta) )  );
                AN = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) ) +  1/(2*sqrt ( ((Uijp-Uimjp)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AC = 1 - AW - AE - AS - AN;


                A( getInd(i,j,m) , getInd(i,j,m) ) = AC;
                A(getInd(i,j,m), getInd(i+1,j,m)) = AE;
                A(getInd(i,j,m), getInd(i,j-1,m)) = AS;
                A (getInd(i,j,m), getInd(i,j+1,m)) = AN;

            elseif i == m 
                Uimj=u(getInd(i-1,j,m));
                Uijm=u(getInd(i,j-1,m));
                Uijp= u (getInd(i,(j+1),m));
                Uipj = 0;
                Uipjm = 0;
                Uimjp = u(getInd(i-1,j+1,m));
                AW = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uimjp-Uimj)/h)^2+beta) )  );
                AE = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uipj-Uipjm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AS = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipjm-Uijm)/h)^2 +((Uij-Uijm)/h)^2+beta) )  );
                AN = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) ) +  1/(2*sqrt ( ((Uijp-Uimjp)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AC = 1 - AW - AE - AS - AN;


                A( getInd(i,j,m) , getInd(i,j,m) ) = AC;
                A( getInd(i,j,m) , getInd(i-1,j,m) ) = AW;
                A(getInd(i,j,m), getInd(i,j-1,m)) = AS;
                A (getInd(i,j,m), getInd(i,j+1,m)) = AN;

            elseif j == 1 
                Uimj=u(getInd(i-1,j,m));
                Uijm=0;
                Uijp= u (getInd(i,(j+1),m));
                Uipj = u(getInd(i+1,j,m));
                Uipjm =0;
                Uimjp = u(getInd(i-1,j+1,m));
                AW = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uimjp-Uimj)/h)^2+beta) )  );
                AE = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uipj-Uipjm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AS = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipjm-Uijm)/h)^2 +((Uij-Uijm)/h)^2+beta) )  );
                AN = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) ) +  1/(2*sqrt ( ((Uijp-Uimjp)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AC = 1 - AW - AE - AS - AN;


                A( getInd(i,j,m) , getInd(i,j,m) ) = AC;
                A( getInd(i,j,m) , getInd(i-1,j,m) ) = AW;
                A(getInd(i,j,m), getInd(i+1,j,m)) = AE;
                A (getInd(i,j,m), getInd(i,j+1,m)) = AN;

            elseif j==m 
                Uimj=u(getInd(i-1,j,m));
                Uijm=u(getInd(i,j-1,m));
                Uijp= 0;
                Uipj = u(getInd(i+1,j,m));
                Uipjm =u(getInd(i+1,j-1,m));
                Uimjp = 0;
                AW = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uimjp-Uimj)/h)^2+beta) )  );
                AE = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uipj-Uipjm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AS = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipjm-Uijm)/h)^2 +((Uij-Uijm)/h)^2+beta) )  );
                AN = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) ) +  1/(2*sqrt ( ((Uijp-Uimjp)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AC = 1 - AW - AE - AS - AN;


                A( getInd(i,j,m) , getInd(i,j,m) ) = AC;
                A( getInd(i,j,m) , getInd(i-1,j,m) ) = AW;
                A(getInd(i,j,m), getInd(i+1,j,m)) = AE;
                A(getInd(i,j,m), getInd(i,j-1,m)) = AS;



            else 
                Uimj=u(getInd(i-1,j,m));
                Uijm=u(getInd(i,j-1,m));
                Uijp= u (getInd(i,(j+1),m));
                Uipj = u(getInd(i+1,j,m));
                Uipjm =u(getInd(i+1,j-1,m));
                Uimjp = u(getInd(i-1,j+1,m));
                AW = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uimjp-Uimj)/h)^2+beta) )  );
                AE = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uipj-Uipjm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AS = -alpha/h^2* (  1/(2*sqrt ( ((Uij-Uimj)/h)^2 +((Uij-Uijm)/h)^2+beta) ) +  1/(2*sqrt ( ((Uipjm-Uijm)/h)^2 +((Uij-Uijm)/h)^2+beta) )  );
                AN = -alpha/h^2* (  1/(2*sqrt ( ((Uipj-Uij)/h)^2 +((Uijp-Uij)/h)^2+beta) ) +  1/(2*sqrt ( ((Uijp-Uimjp)/h)^2 +((Uijp-Uij)/h)^2+beta) )  );
                AC = 1 - AW - AE - AS - AN;


                A( getInd(i,j,m) , getInd(i,j,m) ) = AC;
                A( getInd(i,j,m) , getInd(i-1,j,m) ) = AW;
                A(getInd(i,j,m), getInd(i+1,j,m)) = AE;
                A(getInd(i,j,m), getInd(i,j-1,m)) = AS;
                A (getInd(i,j,m), getInd(i,j+1,m)) = AN;



            end
            
        end
    end 

    mat = A;
end
function ind = getInd(i,j, m)
    ind = i + (j-1) * m ;
end
