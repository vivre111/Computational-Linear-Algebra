time = zeros(4,3);
m = [8,16,24];

for i = 1:3
    [A,b] = Lap2D(m(i));
    tic
    GaussElim(A,b);
    time(i,1) = toc;
    tic 
    Cholesky(A,b);
    time(i,2) = toc;
    BandGE(A,b,m,m);
    time(i,3) = toc;
end
time(4,1)=1000;


m=32;
[A,b] = Lap2D(m);
tic
Cholesky(A,b);
time(4,2)=toc;

tic
BandGE(A,b,m,m);
time(4,3)=toc;

T = array2table(time,'VariableNames',{'GE','Cholesky','BandGE'},'RowNames',{'8','16','24','32'})


