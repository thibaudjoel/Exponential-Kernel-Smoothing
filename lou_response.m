function response_vector = lou_response(h,X,Y)
response_vector = zeros(length(X),1);
kern = @(u) exp(-abs(u));
weights = @(x,X,h) kern(1/h*(X-x));
response_est = @(x,X,Y,h) sum(Y.*weights(x,X,h))/sum(weights(x,X,h));

for i=1:length(X)
    lou_X = X([1:i-1 i+1:end]);
    lou_Y = Y([1:i-1 i+1:end]);
    response_vector(i) = response_est(X(i),lou_X,lou_Y,h);
end
