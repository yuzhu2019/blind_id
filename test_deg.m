clear;clc;
addpath('graph_generator','misc');

N = 30;
% small-world - mean node degree is 4
[A_sw, ~] = generate_connected_smallworld(N,2,0.2);
deg_sw = sum(sum(A_sw));

% ba
temp = [];
for cnt = 1:3000
    seed = generate_connected_ER(6,0.8);
    [A_ba, ~] = generate_connected_BA(N, 2, seed);
    temp = [temp, sum(sum(A_ba))];
end
deg_ba = mean(temp);

% sbm
temp = [];
for cnt = 1:3000
    [A_sbm, ~] = generate_connected_SBM(N,3,0.1,0.21);
    temp = [temp, sum(sum(A_sbm))];
end
deg_sbm = mean(temp);

% ER 
temp = [];
for cnt = 1:3000
    %[A_er, ~] = generate_connected_ER(N,4/30);
    [A_er, ~] = generate_connected_ER(N,0.134);
    temp = [temp, sum(sum(A_er))];
end
deg_er = mean(temp);




