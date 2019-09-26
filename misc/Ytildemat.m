% Y - tilde Y
% N - graph size
% P - number of filters
function Y = Ytildemat(y, V)

    [N, P] = size(y);
    Y = zeros(N*P*(P-1)/2, N*P);
    rows = 1:N;
    colsa = 1:N;
    for i = 1:P
        colsb = (i*N + 1):((i + 1)*N);
        for j = i+1:P
            Y(rows, colsa) =  diag(V'*y(:,j));
            Y(rows, colsb) = -diag(V'*y(:,i));
            rows = rows + N;
            colsb = colsb + N;
        end
        colsa = colsa + N;
    end

end

