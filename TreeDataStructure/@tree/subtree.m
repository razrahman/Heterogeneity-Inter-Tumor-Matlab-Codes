function subtree = subtree(obj, node)
%% SUBTREE  Return the sub-tree made of all the nodes below the given one

    % Get indices of the subtree.
    iterator = obj.depthfirstiterator(node);
    
    % Copy the content of the tree related to the subtree
    parents = obj.Parent(iterator);
    nodes = obj.Node(iterator);
    
    % Revamp parent indices
    newParents = NaN(numel(parents), 1);
    newParents(1) = 0; % The new root node
    for i = 2 : numel(parents)
        pr = parents(i);
        newParents(i) = find( iterator == pr, 1, 'first');
    end
    
    % Create a new tree with the sub-content
    subtree = tree;
    subtree.Node = nodes;
    subtree.Parent = newParents;

end