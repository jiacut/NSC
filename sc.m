function [W, H, cluster_labels] = sc(data, sigma, num_clusters)

% ���������Ծ���W�ͶȾ���D
W = euclidean(data', data');
W = double(exp(-(W.*W) ./ (2*sigma*sigma)));
%  n = size(data,1);	% Ԫ�صĸ���
% W = zeros(n);	% �˾���n*n
% for i = 1:n
%     for j = 1:n
%         dist = sqrt(sum((data(i,:)-data(j,:)).^2)); 
%         W(i,j) = exp(-dist/(2*sigma^2));
%     end
% end
D = sum(W, 2);
D_sqrt = diag(sqrt(1./D)); % D^(-1/2)

% ����Laplacian���������ֽ�
L = D_sqrt * W * D_sqrt;
OPTS.disp = 0;
[V, val] = eigs(L, num_clusters, 'LM', OPTS);

% k-means
sq_sum = sqrt(sum(V.*V, 2)) + 1e-20;
Vn = V ./ repmat(sq_sum, 1, num_clusters);
cluster_labels = kmeans(Vn, num_clusters);

% ������������H
n = size(cluster_labels, 1);
H = zeros(n, num_clusters);
for i = 1:n
    H(i, cluster_labels(i)) = 1;
end


