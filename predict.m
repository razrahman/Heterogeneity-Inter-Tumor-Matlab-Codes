function t = predict(x,t)

if t.isleaf(1)==0
    condition = t.get(1);
    children = t.getchildren(1);
    if x(condition{1})<condition{2}
        subtree = t.subtree(children(1)); % go to left child
    else
        subtree = t.subtree(children(2)); % go to right child
    end
    t = predict(x,subtree);
    
else
    t = t.get(1); % reached the leaf, obtain the avrage output

end