function feature_threshold_t = build_single_tree(X, Y, V_inv, mtree,N,copula,Command, min_leaf, Alpha)

feature_threshold_t=tree;

index=1:size(X,1); % initial index corresponds to all the samples
feature_threshold_t = split_node(X, Y, feature_threshold_t, 1, index, index, V_inv, mtree,N,copula,Command, min_leaf, Alpha);





