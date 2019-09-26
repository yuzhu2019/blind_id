% multi-process, known filter order, noise exits
% influence of the number and lengths of filters
% SBM(30,2,0.1,0.3)
clear;clc;close all;
addpath('graph_generator','misc');

%% parameters
N = 30;                  % graph size
M_vec = [3 3 3 4 6];     % number of filters
L_vec = [5 4 3 3 3];     % filter length
NN = 10.^(-4:0.5:0);     % noise level
num_mc = 2000;           % number of trials    

assert(length(M_vec)==length(L_vec))
%% influence of the number and lengths of filters
err = zeros(num_mc, length(M_vec), length(NN));

for i = 1:num_mc
    disp(i)
    % graph
    [A, ~] = generate_connected_SBM(N,2,0.1,0.3);
    A = A/5;
    [V, Lambda] = eig(A);
    Psi = fliplr(vander(diag(Lambda)));
    % input
    x = randn(N,1);
    
    for j = 1:length(M_vec)
        M = M_vec(j);
        L = L_vec(j);
        % filter
        h = randn(L, M);
        h_true = h(:);
        h_true = h_true/norm(h_true);
        % output without noise
        y_0 = zeros(N, M);
        for m = 1:M
            y_0(:,m) = Hfilt(A, h(:,m))*x;
        end
        % Vandermonde matrix
        SuperPsi = kron(eye(M), Psi(:,1:L));
        
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
            [~, O_S, O_V] = svd(O);
            h_est = O_V(:, end);
            % performance
            h_true = sign(h_true'*h_est)*h_true;
            err(i, j, ni) = norm(h_est-h_true);
        end
    end  
end

%% figure
figure; 
hold on;
plot(NN,mean(squeeze(err(:, 1, :))),'s-','LineWidth',1.5,'MarkerSize',8);
plot(NN,mean(squeeze(err(:, 2, :))),'^-','LineWidth',1.5,'MarkerSize',8);
plot(NN,mean(squeeze(err(:, 3, :))),'o-','LineWidth',1.5,'MarkerSize',8);
plot(NN,mean(squeeze(err(:, 4, :))),'v-','LineWidth',1.5,'MarkerSize',8);
plot(NN,mean(squeeze(err(:, 5, :))),'*-','LineWidth',1.5,'MarkerSize',8);
grid on;
box on;
set(gca, 'Xscale', 'log');
set(gca, 'Yscale', 'log');

xlabel('\sigma_n');
ylabel('Average recovery error');

legend('M=3, L=5', 'M=3, L=4', 'M=3, L=3', 'M=4, L=3', 'M=6, L=3',...
       'Location', 'NorthWest');

hold off;
