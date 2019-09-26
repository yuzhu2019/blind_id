% multi-process, unknown filter order, noise exits
% influence of overshot filter lengths and weight matrices
clear;clc;close all;
addpath('graph_generator','misc');

%% parameters
N = 30;
p = 0.1;
M = 3;
L = 3;
Q = 5; 
NN = 10.^(-4:0.5:0);
num_mc = 1000;
% W = eye(M*Q-1);
W = exp_wmat(M, Q, 1); 

%% algorithms
err = zeros(length(NN), num_mc);

for i = 1:num_mc
    disp(i)
    % graph  
    [A, ~] = generate_connected_ER(N, p);
    [V, Lambda] = eig(A);
    Psi = fliplr(vander(diag(Lambda)));
    SuperPsi = kron(eye(M), Psi(:,1:Q));
    % input
    x = randn(N, 1);
    % filter
    h = randn(L, M);
    h(1,1) = 1;
    hbar = [h; zeros(Q-L, M)];
    hbar = hbar(:);
    h_true = hbar(2:end);
    h_true_2 = norm(h_true);
    % output without noise
    y_0 = zeros(N, M);
    for m = 1:M
        y_0(:,m) = Hfilt(A, h(:,m))*x;
    end
    
    for ni = 1:length(NN)
        % add noise        
        y = zeros(N, M);
        for m = 1:M
            sigma = sqrt(norm(y_0(:,m))^2/N) * NN(ni);
            y(:,m) = y_0(:,m) + sigma * randn(N, 1); 
        end
        % algorithm
        Y = Ytildemat(y, V);
        O = Y*SuperPsi;
        b = O(:, 1);
        Phi = O(:, 2:end);
        cvx_begin quiet
            variable h_est(M*Q-1, 1)
            minimize norm(W*h_est, 1)
            subject to
                norm(Phi*h_est+b) <= norm(Phi*h_true+b);
        cvx_end
        % performance
        err(ni,i) = norm(h_est-h_true)/h_true_2;
    end
end


