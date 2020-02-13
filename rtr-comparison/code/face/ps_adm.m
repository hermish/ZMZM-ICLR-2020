function p = ps_adm(Y, p_init, tau, tol, max_iter) 

done = false;
iter = 0; 

p_old = p_init; 

while ~done,    
	iter = iter + 1;          
    x = soft_thresholding(p_old' * Y, tau);    
    p_new = proj_sphere(Y * x');       
    
    if norm(p_old - p_new) < tol || iter > max_iter
    	done = true; 
    end 
    
    disp([' Iter: ' num2str(iter) ' Obj: ' num2str(tau*norm1(x) + 0.5*norm(p_new' * Y - x)^2) ' Movement: ' num2str(norm(p_old - p_new)) ])
    
    p_old = p_new; 
end

p  = p_new; 
