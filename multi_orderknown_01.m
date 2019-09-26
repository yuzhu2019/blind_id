% multi-process, known filter order, noiseless
% verification of Theorem 1 and 2
clear;clc;close all;
addpath('graph_generator','misc');

%% parameters
N = 25;                     % graph size
p = 0.2;                    % edge prob.
L = 8;                      % filter length
M = 5;                      % number of filters
OMEGA = 1:24;               % input richness
a = [0, 0.5, 0.8, 0.9];     % filter correlation
num_mc = 1000;              % number of trials

%% verify Theorem 1 and 2
err = zeros(num_mc, length(a), length(OMEGA));

for i = 1:num_mc
    disp(i)
    % graph  
    [A, ~] = generate_connected_ER(N, p);
    [V, Lambda] = eig(A);
    % Vandermonde matrix
    Psi = fliplr(vander(diag(Lambda)));
    SuperPsi = kron(eye(M), Psi(:,1:L));
    
    for j = 1:length(a)
        % filter
        h = randn(L, M);
        h(:, 2:end) = a(j)*repmat(h(:,1), 1, M-1) + (1-a(j))*randn(L, M-1);
        % filter matrix
        H = zeros(N, N, M);
        for m = 1:M
            H(:,:,m) = Hfilt(A, h(:,m));
        end
        h_true = h(:);
        h_true = h_true/norm(h_true);
        
        for o = 1:length(OMEGA)
            % input
            omega = OMEGA(o);
            xtilde = zeros(N, 1);
            xtilde(1:omega) = randn(omega, 1);
            x = V*xtilde;  
            % output
            y = zeros(N, M);
            for m = 1:M
                y(:,m) = H(:,:,m)*x;
            end
            % algorithm
            Y = Ytildemat(y, V);
            O = Y*SuperPsi;
            [~, O_S, O_V] = svd(O);
            h_est = O_V(:, end);
            % performance
            h_true = sign(h_true'*h_est)*h_true;
            err(i,j,o) = norm(h_est-h_true);
        end
    end
end

%% figure
mean_err = squeeze(mean(err < 0.01));
figure; hold on;
plot(OMEGA, mean_err(1,:), 's-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(OMEGA, mean_err(2,:), '^-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(OMEGA, mean_err(3,:), 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(OMEGA, mean_err(4,:), 'x-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot([9 9], [-0.1 1.1], '--', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2]);
plot([15 15], [-0.1 1.1], '--', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2]);

axis([OMEGA(1) OMEGA(end) -0.1 1.1]);           

xlabel('|\Omega_1|');
ylabel('Recovery rate');

grid on;
box on;

legend('a=0', 'a=0.5', 'a=0.8', 'a=0.9', ...
       'Location', 'SouthEast');

