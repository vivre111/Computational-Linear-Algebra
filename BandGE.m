function x = BandGE(A,b,plow,qup)
s = length(A);
for j = 1:(s-1)
    for i = min(j+qup,s):-1:j+1
        m = A(i,j)/A(j,j);
        A(i,max(j-plow,1):min(j+qup,s)) = A(i,max(j-plow,1):min(j+qup,s)) - m*A(j,max(j-plow,1):min(j+qup,s));
        b(i) = b(i) - m*b(j);
    end
end 
x = zeros(s,1);
x(s) = b(s)/A(s,s);               
for i = s-1:-1:1                    
    x(i) = (b(i)- sum(A(i,i+1:min(i+qup,s))*x(i+1:min(i+qup,s))) )/A(i,i);

end 
end
