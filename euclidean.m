function d = euclidean(a, b)
%EUCLIDEAN Compute the Euclidean distance matrix between two matrices.
%   This function is designed for processing very large data using divide-
%   and-conquer technique.
%
%   Input : a : D-by-M data matrix, where D is the number of dimensions,
%               M is the number of data
%           b : D-by-N data matrix, where D is the number of dimensions,
%               N is the number of data
%   Output: d : M-by-N matrix of the Euclidean distance between a and b.
%欧氏距离矩阵计算两个矩阵之间的欧氏距离矩阵。
%此函数用于使用divide处理非常大的数据-
%征服技术。
%
%输入：a:D-by-M数据矩阵，其中D是维度数，
%        M是数据的数量
%       b:D-x-N数据矩阵，其中D是维数，
%        N是数据的数量
%输出：d: a和b之间的欧几里得距离的M×N矩阵。
%
% Calculate a^2, b^2, here we assume b is larger than a：这里我们假设b比a大
%

aa = single(full(sum(a.*a, 1)));
bb = single(full(sum(b.*b, 1)));

%
% Do a*b in several steps instead of once because of memory limitation，由于内存限制，a*b分几个步骤执行，而不是一次
%
two_ab = single(zeros(size(aa, 2), size(bb, 2)));
% Select at most 10000 instances of b for a*b per iteration
num_iter = ceil(size(bb, 2)/10000);
for i = 1:num_iter
  start_index = 1 + (i-1)*10000;
  end_index = min(i*10000, size(bb, 2));
  abtmp = single(full(a'*b(:, start_index:end_index)));
  two_ab(:, start_index:end_index) = 2*abtmp;
end % Now we have entire ab
clear a b abtmp;

d = bb(ones(size(aa, 2), 1), :);
d = d - two_ab; % Now we have d = b^2 - 2ab
clear two_ab;

ff = aa';
ff = ff(:, ones(size(bb, 2), 1));
d = d + ff; % Now we have d = a^2 + b^2 -2ab
clear aa bb ff;
d = sqrt(d);