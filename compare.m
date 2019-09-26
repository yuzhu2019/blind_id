% single, unknown, noise exists
clear;clc; 
close all;
addpath('misc');
load('karate.mat');

%% graph
N = 34;
G = graph(edges(:,1), edges(:,2));
A = full(adjacency(G));  
A = A/3;
[V, Lambda] = eig(A); % non-zero: 24
Psi = fliplr(vander(diag(Lambda)));

%% parameters
L_vec = [3 5 7];
Q_vec = [7 7 7];

assert(length(L_vec)==length(Q_vec))
M = length(L_vec);  % number of filters
L_max = max(L_vec); % filter total length
Q_max = max(Q_vec);

NN = 10.^(-6:1:0);
num_mc = 1000;
%% leverage decreasing coefficients info (single, unknown)
temp = eye(Q_max);
P2 = [];
for m = 1:M
    P2 = [P2, temp(:, 1:Q_vec(m))];
end
P1 = zeros(Q_max-1, Q_max);
for m = 1:Q_max-1
    P1(m,m) = 1;
    P1(m,m+1) = -1;
end
P = P1*P2;

%% multi-process, order known/unknown
SuperPsi = [];
SuperTheta = [];
for m = 1:M
    SuperPsi = blkdiag(SuperPsi, Psi(:,1:L_vec(m)));
    SuperTheta = blkdiag(SuperTheta, Psi(:,1:Q_vec(m)));
end

%% single process, order known/unknown
SuperPsiBar = zeros(N*M, L_max);
SuperThetaBar = zeros(N*M, sum(Q_vec));
temp = [];
rows = 1:N;
for m = 1:M
    temp = cat(2,temp,Psi(:,1:Q_vec(m)));
    SuperPsiBar(rows, 1:L_vec(m)) = Psi(:, 1:L_vec(m));
    SuperThetaBar(rows, 1:sum(Q_vec(1:m))) = temp;
    rows = rows + N;
end

%% algorithm
h_l = 0.2;
h_u = 1;

err1 = zeros(num_mc, length(NN));
err2 = zeros(num_mc, length(NN));
err3 = zeros(num_mc, length(NN));
err4 = zeros(num_mc, length(NN));
err5 = zeros(num_mc, length(NN));
err6 = zeros(num_mc, length(NN));

for i = 1:num_mc
    disp(i)
    % input
    x = randn(N,1);
    % filter
    h = h_l + (h_u - h_l)*rand(L_max,1); % uniformly distributed (h_l, h_u)
    h = sort(h, 'descend');    
    h(1) = 1;
    % output without noise
    y_0 = zeros(N,M);
    for m = 1:M
        y_0(:,m) = Hfilt(A, h(1:L_vec(m)))*x;
    end
    
    %% true filter coefficients
    h_true1 = [];
    h_true2 = [];
    for m = 1:M
        h_true1 = cat(1,h_true1,h(1:L_vec(m)));
        h_true2 = cat(1,h_true2,[h(1:L_vec(m));zeros(Q_vec(m)-L_vec(m),1)]);
    end
    h_true4 = [h(1:L_vec(1)); zeros(Q_vec(1)-L_vec(1),1)];
    for m = 2:M
        temp = [zeros(L_vec(m-1),1); h(L_vec(m-1)+1:L_vec(m),1); zeros(Q_vec(m)-L_vec(m),1)];
        h_true4 = cat(1,h_true4,temp);
    end
  
    h_true1 = h_true1/norm(h_true1);
    h_true2 = h_true2(2:end);
    h_true3 = h/norm(h);
    h_true4 = h_true4(2:end);
   
    %% different noise levels
    for ni = 1:length(NN)
        % add noise
        y = zeros(N,M);
        for m = 1:M
            sigma = sqrt(norm(y_0(:,m))^2/N) * NN(ni);
            y(:,m) = y_0(:,m) + sigma * randn(N, 1); 
        end
        Y = Ytildemat(y, V);
        
        %% (1) algorithm - multi-process, order known
        O1 = Y*SuperPsi;
        [~, O1_S, O1_V] = svd(O1);
        h_est1 = O1_V(:, end);
        
        h_true1 = sign(h_true1'*h_est1)*h_true1;
        err1(i,ni) = norm(h_est1-h_true1);
        
        %% (2) algorithm - multi-process, order unknown
        O2 = Y*SuperTheta;
        b = O2(:, 1);
        Phi = O2(:, 2:end);
        cvx_begin quiet
            variable h_est2(sum(Q_vec)-1, 1)
            minimize norm(h_est2, 1)
            subject to
                norm(Phi*h_est2+b) <= norm(Phi*h_true2+b);
        cvx_end
      
        err2(i,ni) = norm(h_est2-h_true2)/norm(h_true2);
        
        %% (3) algorithm - single process, order known
        O3 = Y*SuperPsiBar;
        [~, O3_S, O3_V] = svd(O3);
        h_est3 = O3_V(:, end);
        
        h_true3 = sign(h_true3'*h_est3)*h_true3;
        err3(i,ni) = norm(h_est3-h_true3);
        
        %% (4) algorithm - single process, order unknown
        O4 = Y*SuperThetaBar;
        b = O4(:, 1);
        Phi = O4(:, 2:end);
        cvx_begin quiet
            variable h_est4(sum(Q_vec)-1, 1)
            minimize norm(h_est4, 1)
            subject to
                norm(Phi*h_est4+b) <= norm(Phi*h_true4+b);
        cvx_end
        
        err4(i,ni) = norm(h_est4-h_true4)/norm(h_true4);
       
        
        %% with constraint 1
        cvx_begin quiet
            variable h_est4(sum(Q_vec)-1, 1)
            minimize norm(h_est4, 1)
            subject to
                norm(Phi*h_est4+b) <= norm(Phi*h_true4+b);
                h_est4 >= 0;
        cvx_end
        
        err5(i,ni) = norm(h_est4-h_true4)/norm(h_true4);
        
        %% with constraint 2
        cvx_begin quiet
            variable h_est4(sum(Q_vec)-1, 1)
            minimize norm(h_est4, 1)
            subject to
                norm(Phi*h_est4+b) <= norm(Phi*h_true4+b);
                P(:,2:end) * h_est4 + P(:,1) >= 0;
                h_est4 >= 0;
        cvx_end
        
        err6(i,ni) = norm(h_est4-h_true4)/norm(h_true4);
        
    end
end




