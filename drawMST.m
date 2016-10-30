function [ ] = drawMST(treeNodes, treeEdges, nodes, nodeStatus)
% draw nodes matrix in figure
scatter(treeNodes(1,2),treeNodes(1,1), 'Marker', 'square', 'MarkerEdgeColor', 'red'); 
axis ij; axis([1, nodes(2), 1, nodes(1)]); hold on;
scatter(treeNodes(2:end,2),treeNodes(2:end,1), 'MarkerEdgeColor', 'blue'); hold off;
%% draw edge for each node pair in tree
classes = unique(nodeStatus);
colorList = {'cyan', 'red', 'yellow', 'green', 'blue', 'magenta', 'black', ...
    [0.6 0.2 0], [0.5 0.5 0.5]};
for i = 1 : size(treeEdges, 1)
   x1 = treeEdges(i, 2);
   x2 = treeEdges(i, 4);
   y1 = treeEdges(i, 1);
   y2 = treeEdges(i, 3);
   index = find(classes == nodeStatus(y1, x1));
   line([x1 x2], [y1 y2], 'Color', colorList{index});
   text((x1+x2)/2, (y1+y2)/2, num2str(treeEdges(i,5)), 'Color', colorList{index});
end

end

