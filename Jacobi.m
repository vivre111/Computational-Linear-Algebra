function [x,iter] = Jacobi(A,b,x_initial,maxiter,tol)
    
    x=x_initial;
    iter = 0;
    for i = 1:maxiter
        r = b - A*x;
        if norm(r) < tol * norm(b)
            return
        end
        iter = iter + 1;
        x = x + diag( diag(A))^(-1) * r;
    end

end 

