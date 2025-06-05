% Function to create Regression Matrix
function X = create_regression_matrix(y, u, na, nb, nk)
        % y: Output data
        % u: Input data
        % na: Order of A(q)
        % nb: Order of B(q)
        % nk: Dead time (delay)
 N = length(y); % No. of Samples
 X = zeros(N - na, na + nb); % Initialize Regression matrix

% Construct ARX regression matrix
    for i = 1:N
        % fill in past output values (AR part)
        for j = 1:na
            if (i - j) > 0
                X(i, j) = -y(i - j);
            end
        end
        % fill in past input values (X part)
        for j = 1:nb
            if(i -nk -j +1) > 0
                X(i, na + j) = u(i -nk - j + 1);
            end
        end
    end
end
