% generate exponential weight matrix
function W = exp_wmat(M, Q, alpha)

w = repmat(exp(alpha*(1:Q)), 1, M);
w = w(2:end); % remove the first entry
W = diag(w);

end