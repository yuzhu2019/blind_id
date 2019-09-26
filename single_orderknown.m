% single process, known filter order, noise exists
clear;clc;
close all;
addpath('graph_generator','misc');

%% parameters
N = 30;                
p = 0.1;               
w_l = 0.1;             
w_u = 0.7;             
L_max = 8;
LLc = {[2 4 6 8], [3 5 8], [4 8]};
NN = 10.^(-5:0.5:0);
num_mc = 2000;

%% algorithms
err = zeros(length(LLc), length(NN), num_mc);

for i = 1:num_mc
    disp(i)
    % graph
    [A, ~] = generate_connected_ER(N, p); 
    [A, ~] = add_weights(A, w_l, w_u);
    [V, Lambda] = eig(A);
    Psi = fliplr(vander(diag(Lambda)));
    % input
    x = randn(N, 1);  
    % filter
    h = randn(L_max, 1);   
    h_true = h/norm(h);
    % output without noise
    y_0 = zeros(N, L_max);
    for m = 1:L_max
        y_0(:,m) = Hfilt(A, h(1:m)) * x;
    end
    
    for ni = 1:length(NN)
        % add noise
        y_all = zeros(N, L_max);
        for m = 1:L_max
            sigma = sqrt(norm(y_0(:,m))^2/N) * NN(ni);
            y_all(:,m) = y_0(:,m) + sigma * randn(N, 1); 
        end
        
        for j = 1:length(LLc)
            % algorithm
            L = LLc{j};
            M = length(L);
            y = y_all(:,L);
            Y = Ytildemat(y, V);
            SuperPsiBar = zeros(M*N, L(end));
            rows = 1:N;
            for m = 1:M
                cols = 1:L(m);
                SuperPsiBar(rows,cols) = Psi(:,cols);
                rows = rows + N;
            end
            O = Y*SuperPsiBar;
            [~, O_S, O_V] = svd(O);
            h_est = O_V(:, end);
            % performance
            h_true = sign(h_true'*h_est)*h_true;
            err(j,ni,i) = norm(h_est-h_true);
        end
    end
end

%% figure
figure;
hold on;
plot(NN, mean(squeeze(err(3,:,:)),2), 's-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(NN, mean(squeeze(err(2,:,:)),2), '^-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(NN, mean(squeeze(err(1,:,:)),2), 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);

grid on;
box on;

set(gca, 'Xscale', 'log');
set(gca, 'Yscale', 'log');

xlabel('\sigma_n');
ylabel('Average recovery error');

legend('L = [4 8]', 'L = [3 5 8]', 'L = [2 4 6 8]', ...
       'Location', 'SouthEast');

