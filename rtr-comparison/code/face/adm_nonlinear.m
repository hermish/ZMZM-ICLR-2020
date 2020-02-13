function [q,y] = adm_nonlinear(A,q_init,lambda,MaxIter,tol)

q = q_init;
for k = 1:MaxIter
    q_old = q;
    y = soft_thresholding(A*q,lambda);
    q = A'*y/norm(A'*y,2);
    res_q = norm(q_old-q,2);
    if (res_q<=tol)
        return;
    end 
    
%     if(mod(k,10)==0)
%         fprintf('Running the %d-th iteration, diff_q=%f, diff_obj =%f \n',k, res_q,res_obj);
%     end
end

end

function Y = soft_thresholding(X,d)

Y = sign(X).*max(abs(X)-d,0); 
end

% function Y = soft(X,rho)
%     Y = (X-rho*sign(X)).*(abs(X)>=rho);
% end