function theta = estimate_arx(X,y)
         % X: Regression matrix
         % y: Output data
    theta = (X' * X) \ (X' * y);
end