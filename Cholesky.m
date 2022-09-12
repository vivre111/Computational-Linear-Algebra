function x = Cholesky(A,b)
    G = choleskyMat(A);
    temp = forward(G,b);
    x = backward(transpose(G), temp);
end

function [F]=choleskyMat(A)
    [m,n]=size(A);
    L=zeros(m,m);%Initialize to all zeros
    row=1;col=1;
    j=1;
    for i=1:m
        a11=sqrt(A(1,1));
        L(row,col)=a11;
        if(m~=1) %Reached the last partition
            L21=A(j+1:m,1)/a11;
            L(row+1:end,col)=L21;
            A=(A(j+1:m,j+1:m)-L21*L21');
            [m,n]=size(A);
            row=row+1;
            col=col+1;
        end
    end
    F=L;
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


function[x]=backward(U,b)
S=size(U);
m=S(1);
x=zeros(1,m);
x(1,m)=b(end)./U(m,m);
%bacward substitution
for k=m-1:-1:1
        x1=1/U(k,k).*(b(k)-sum(U(k,k+1:end).*x(k+1:end)));
        x(k)=x1;
end
x=x';
end
