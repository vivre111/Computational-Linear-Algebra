function [x,iter] = SOR(omega, A,b,x_initial,maxiter,tol)
    
    x=x_initial;
    iter = 0;
    for i = 1:maxiter
        r = b - A*x;
        if norm(r) < tol * norm(b)
            return
        end
        iter = iter + 1;
        L =  ( 1/omega * diag( diag(A))- tril(A,-1));
        % we want x=L^(-1) r, instead we solve x s.t. Lx = r
        x = x + forward(L,r);
    end

end 



function[x]=forward(L,b)
S=size(L);
m=S(1);
x=zeros(1,m);
x(1,1)=b(1)./L(1,1);
for k=2:m  
        x1=1/L(k,k).*(b(k)-sum(L(k,k-1:-1:1).*x(k-1:-1:1)));
        x(1,k)=x1;
end
x=x';
end
