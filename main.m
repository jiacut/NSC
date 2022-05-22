clear
% clc

[data, label] = load_dataset(1); % 加载数据

num_points = length(label);
num_clusters = length(unique(label));
sigma =0.5;

alpha =0.26; 
beta =10;

% Compute the kernel
    [K, H, cluster_labels_sc] = sc(data, sigma, num_clusters);
    result_sc = ClusteringMeasure(label, cluster_labels_sc);


% 欧式距离的映射
    one_n = ones(num_points, num_points);
    E = one_n * diag(diag(K)) + diag(diag(K)) * one_n - 2 * K;

% 非负子空间核聚类
    [cluster_labels, epoch] = nskc(alpha, beta, K, H, E);
    result = ClusteringMeasure(label, cluster_labels);





