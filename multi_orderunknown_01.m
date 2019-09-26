% multi-process, unknown filter order, noiseless
% verification of Theorem 4
clear;clc;close all;
addpath('graph_generator','misc');

%% parameters
N = 30;                        % graph size
p = 0.1;                       % edge prob.
M = 3;                         % number of filters
L = 3;                         % filter length
Q_vec = [4 6 6];               % supposed filter length
W_vec = [false, false, true];  % false-identify, true-exp weighting
num_mc = 20000;
delta = 0.02;

assert(length(Q_vec)==length(W_vec)) % number of settings
p_options = optimoptions('linprog', 'Display', 'off');

%% verification of Theorem 4
cond1 = false(length(Q_vec),num_mc);
cond2 = zeros(length(Q_vec),num_mc);
err = zeros(length(Q_vec),num_mc);

for i = 1:num_mc
    disp(i)
    % graph  
    [A, ~] = generate_connected_ER(N, p);
    [V, Lambda] = eig(A);
    Psi = fliplr(vander(diag(Lambda)));
    % input
    x = randn(N, 1);
    % filter
    h = randn(L, M);
    h(1,1) = 1;
    % output
    y = zeros(N, M);
    for m = 1:M
        y(:,m) = Hfilt(A, h(:,m))*x;
    end
    Y = Ytildemat(y, V);
    
    for j = 1:length(Q_vec)
        % Phi & b
        Q = Q_vec(j);
        SuperPsi = kron(eye(M), Psi(:,1:Q));
        O = Y*SuperPsi;
        b = O(:, 1);
        Phi = O(:, 2:end);
        % support set
        hbar = [h; zeros(Q-L, M)];
        hbar = hbar(:);
        h_true = hbar(2:end);
        sI = find(h_true);
        sJ = find(h_true == 0);
        % check condition (i)
        cond1(j,i) = rank(Phi(:,sI)) == length(sI);
        % weight matrix
        W = eye(M*Q-1);
        if W_vec(j) == true
            W = exp_wmat(M, Q, 1);
        end
        W_inv = diag(1./diag(W));
        % condition (ii), compute psi
        I = eye(M*Q-1);
        psimat = I(:,sJ)' * inv(delta^-2*W_inv*Phi'*Phi*W_inv + I(:,sJ)*I(:,sJ)') * I(:,sI);
        cond2(j,i) = norm(psimat,Inf);
        % algorithm        
        p_f = [diag(W); zeros(M*Q-1, 1)];
        p_A = [-I,I;-I,-I];
        p_b = zeros(2*(M*Q-1), 1);
        p_Aeq = [zeros(size(Phi)), Phi];
        p_beq = -b;
        ret = linprog(p_f, p_A, p_b, p_Aeq, p_beq, [], [], p_options);
        h_est = ret((M*Q):end);
        % performance
        err(j,i) = norm(h_est-h_true)/norm(h_true);
    end
end

success = err < 0.01;
%% figure
figure
subplot(311)
hold on
histogram(cond2(1, success(1, :)),  linspace(0, 5, 21))
histogram(cond2(1, ~success(1, :)), linspace(0, 5, 21))
vline(median(cond2(1,:)), 'k--')
vline(median(cond2(1, success(1, :))), 'b--')
vline(median(cond2(1, ~success(1, :))), 'r--')
vline(1,'-r')
hold off
axis([-0.25 5.25 0 1500]);   
legend({'Success', 'Failure'})
xlabel('\psi')
ylabel('Frequency')
title(sprintf('No weighting Q=%d', Q_vec(1)))
fprintf('No weighting Q=%d: success prob=%.2f.\n', Q_vec(1), sum(success(1,:))/num_mc)

subplot(312)
hold on
histogram(cond2(2, success(2, :)),  linspace(0, 5, 21))
histogram(cond2(2, ~success(2, :)), linspace(0, 5, 21))
vline(median(cond2(2,:)), 'k--')
vline(median(cond2(2, success(2, :))), 'b--')
vline(median(cond2(2, ~success(2, :))), 'r--')
vline(1,'-r')
hold off
axis([-0.25 5.25 0 1500]); 
legend({'Success', 'Failure'})
xlabel('\psi')
ylabel('Frequency')
title(sprintf('No weighting Q=%d', Q_vec(2)))
fprintf('No weighting Q=%d: success prob=%.2f\n', Q_vec(2), sum(success(2,:))/num_mc)

subplot(313)
hold on
histogram(cond2(3, success(3, :)),  linspace(0, 5, 21))
histogram(cond2(3, ~success(3, :)), linspace(0, 5, 21))
vline(median(cond2(3,:)), 'k--')
vline(median(cond2(3, success(3, :))), 'b--')
vline(median(cond2(3, ~success(3, :))), 'r--')
vline(1,'-r')
hold off
axis([-0.25 5.25 0 6000]); 
legend({'Success', 'Failure'})
xlabel('\psi')
ylabel('Frequency')
title(sprintf('Exponential weighting Q=%d', Q_vec(3)))
fprintf('Weighting Q=%d: success prob=%.2f\n', Q_vec(3), sum(success(3,:))/num_mc)


