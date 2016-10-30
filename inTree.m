function status = inTree(node, nodeStatus)
    if nodeStatus(node(1), node(2)) == 1  % node status = 1, in tree; = 0, not in tree
        status = true;
    else
        status = false;
    end
end