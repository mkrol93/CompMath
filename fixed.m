function fixed(g,p0,tol,n)
% This function solves p=g(p) undercertain conditions on g.
iter=0;
p1=g(p0);
err=abs(p1-p0);
disp('__________________________________________')
disp('iter         pn         g(pn)          err')
disp('__________________________________________')
fprintf('%.0f  %10.6f  %10.6f  %.6f\n',iter,p0,p1,err)
while (err>tol)&&(iter<=n)
    p0=p1;
    p1=g(p0);
    iter=iter+1;
    err=abs(p1-p0);
    fprintf('%.0f  %10.6f  %10.6f  %.6f\n',iter,p0,p1,err)
end
if iter>n
    disp('The sequence do not converge.')
end
end

