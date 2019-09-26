% Add weights (randomly selected from U(a,b)) to an unweighted graph
% Input: unweighted adjacency matrix A_in, range of weights a,b
% Ouput: weighted adjacency matrix/combinatorial Laplacian A_out/L_out
% Yu Zhu, Rice ECE, 12/04/2018
function [A_out,L_out] = add_weights(A_in,a,b)
    N_v = size(A_in,1);
    weights = a + (b-a)*rand(N_v); 
    A_out = A_in .* weights;
    A_out = triu(A_out,1);
    A_out = A_out + A_out';
    L_out = diag(sum(A_out))-A_out;
end

