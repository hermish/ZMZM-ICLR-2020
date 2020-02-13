function X = soft_thresholding(X,d)

X = sign(X) .* max( abs(X) - d, 0 ); 


