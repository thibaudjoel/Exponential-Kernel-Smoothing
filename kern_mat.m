function kernel_matrix = kern_mat(h,X)
n = length(X);
kernel_matrix = zeros(n,n);
kern = @(u) exp(-abs(u));
weights = @(x,h) kern(1/h*(X-x));
for i=1:n
    kernel_matrix(i,:) = weights(X(i),h);
    kernel_matrix(i,:) = kernel_matrix(i,:)/sum(weights(X(i),h));
end