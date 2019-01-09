function ourjacobi(A,b,x0,tol,itmax)
% This MATLAB function solves Ax=b by Jacobi's method.
n=length(b);
for i=1:n
    if A(i,i)==0
        disp('Jacobi does not work.')
        return
    end
end
D=diag(diag(A));
L=tril(A,-1);
U=triu(A,1);
Tjac=-inv(D)*(L+U);
cjac=inv(D)*b;
x=x0;
err=1000;
k=1;
while k<=itmax && err>tol
    x=Tjac*x+cjac;
    err=norm(x-x0);
    x0=x;
    k=k+1;
end
if k==itmax+1
    disp('No convergence. Increase itmax.')
end
fprintf('\n')
disp(['The approximate solution after ',num2str(k-1),'iterations is'])
x   
end

