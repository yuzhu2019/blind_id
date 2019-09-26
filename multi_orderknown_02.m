% multi-process, known filter order, noise exits
% influence of the input spectral richness
% small world graph (30,2,0.2)
clear;clc;%close all;
addpath('graph_generator','misc');

%% parameters
N = 30;                 % graph size
L = 3;                  % filter length
M = 3;                  % number of filters
OMEGA = 12:6:24;        % input spectral richness
NN = 10.^(-5:0.5:-1);   % noise level
num_mc = 2000;          % number of trials

%% influence of the input spectral richness
err = zeros(num_mc, length(OMEGA), length(NN));

for i = 1:num_mc
    disp(i)
    % graph  
    %[A, ~] = generate_connected_ER(N,4/30);
    %[A, ~] = generate_connected_SBM(N,3,0.1,0.21);
    [A, ~] = generate_connected_smallworld(N,2,0.2);
    %seed = generate_connected_ER(6,0.8);
    %[A, ~] = generate_connected_BA(N, 2, seed);
    
    [V, Lambda] = eig(A);
    % Vandermonde matrix
    Psi = fliplr(vander(diag(Lambda)));
    SuperPsi = kron(eye(M), Psi(:,1:L));
    % filter
    h = randn(L, M);
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
        % output without noise
        y_0 = zeros(N, M);
        for m = 1:M
            y_0(:,m) = H(:,:,m)*x;
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
            [~, O_S, O_V] = svd(O);
            h_est = O_V(:, end);
            % performance
            h_true = sign(h_true'*h_est)*h_true;
            err(i, o, ni) = norm(h_est-h_true);
        end
    end
end

%% figure
%figure; 
hold on;
plot(NN,mean(squeeze(err(:, 1, :))),'s-','LineWidth',1.5,'MarkerSize',8);
plot(NN,mean(squeeze(err(:, 2, :))),'^-','LineWidth',1.5,'MarkerSize',8);
plot(NN,mean(squeeze(err(:, 3, :))),'o-','LineWidth',1.5,'MarkerSize',8);
grid on;
box on;
set(gca, 'Xscale', 'log');
set(gca, 'Yscale', 'log');

xlabel('\sigma_n');
ylabel('Average recovery error');

legend('|\Omega_1|=12', '|\Omega_1|=18', '|\Omega_1|=24', ...
        'Location', 'SouthEast');

%hold off;




