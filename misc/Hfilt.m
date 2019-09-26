% generate the filter matrix H 
% given the GSO S and the filter coefficients h

function H = Hfilt(S, h)
    N = size(S, 1);
    H = zeros(N, N);
    for l = 1:length(h)
        H = H + h(l)*S^(l-1);
    end    
end
