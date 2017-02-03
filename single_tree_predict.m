function Y = single_tree_predict(X_hat, t, variable_num)
% X: input testing sets
% t: tree model built
% Y: output predicted value

Y = zeros(size(X_hat,1),variable_num+1);
for ii = 1:size(X_hat,1)
    x = X_hat(ii,:);
    
    leaf_info = predict(x,t);
    Y(ii,:) = [mean(leaf_info(:,2:end-1),1) mode(leaf_info(:,end))];
end

