function q = admm(A,q_init,MaxIter,tol)

rho = 0.1;
q = q_init;
y = A*q;
tau = y;
obj = inf;
for k = 1:MaxIter
    y_old = y;
    y = soft_thresholding(A*q-tau/rho,1/rho);
    q = A'*(rho*y+tau);
    q = q/norm(q,2);
    tau = tau + rho*(y-A*q);
    res_p = norm(y-A*q);
    res_d = norm(rho*(y-y_old),2);
    if (res_p>=10*res_d)
        rho = rho*5;
    end
    if (res_d>=10*res_p)
        rho = rho/5;
    end
    
    if (res_p<=tol&res_d<=tol)
        return;
    end
    
    if(mod(k,10)==0)
        fprintf('Running the %d-th simulation: Primal Res = %f Dual Res = %f\n',...
            k,res_p,res_d);
    end
end

end

function Y = soft_thresholding(X,d)
Y = sign(X).*max(abs(X)-d,0);
end
