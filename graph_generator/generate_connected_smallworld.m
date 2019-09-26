function [A,L] = generate_connected_smallworld(N,K,beta)

    flag_connected = 0;
    while flag_connected ==0
        h = WattsStrogatz(N,K,beta);
        A = full(adjacency(h));
        L = diag(sum(A)) - A; 
        if sum(abs(eig(L)) <= 1.00e-06) == 1
            flag_connected = 1; % Check whether the graph is connected
        end
    end
    
end
