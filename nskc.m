function [cluster_labels, epoch] = nskc(alpha, beta, K, H, E)

epoch = 0;
max_epoch = 1000;
gap = 1;
min_gap = 1e-5;

S = E;
D_sqrt = diag(sqrt(1 ./ sum(S, 2)));
H = H + 0.5;
while (gap > min_gap) && (epoch < max_epoch)
    % update S
    P = D_sqrt * H * H' * D_sqrt;
    S = S .* sqrt((K + beta * P) ./ (K * S + alpha * E)); % tr(S^T * E)
    
    % update H
    S = (S + S') / 2;
    D_sqrt = diag(sqrt(1 ./ sum(S, 2)));
    Q = D_sqrt * S * D_sqrt * H;
    Ht = H .* sqrt(Q ./ (H * H' * Q));

	gap = sum(sum((Ht- H).^2));
	H = Ht;
    epoch = epoch + 1;
end

[~, cluster_labels] = max(H, [], 2);

