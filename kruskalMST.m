%% thie script generate a minimum spanning tree from a downsampled image
% adjustable parameters:
%     scale = 4, 2, 1
%     neighbor = 8, 4 

clear all;
% close all;
clc;

scale = 1;      % 4, 2 or 1 dowbsample ratio
neighbor = 8;   % 8 or 4 neighbors
fprintf('Set parameters: scale = %d, neighbor = %d\n', scale, neighbor);

imgFile = './5_1.jpg';  % load a gray image
img = imread(imgFile);
if ndims(img) > 2
    img = img(:, :, 1);
end
img = single(img); 
dims = size(img);

figure;
subplot(2, 2, 1);
imagesc(img); title(['original image, ', num2str(dims(1)), ' x ', num2str(dims(2))]); axis off;
fprintf('Using image file: %s\n', imgFile);
%% calculate neighbor costs in MST model
if neighbor == 8    % define the 8-neighbor set
    dn_area = [1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1;  1 -1];
elseif neighbor == 4    % define the 4-neighbor set
    dn_area = [1 0; 0 1; -1 0; 0 -1];
end
nodes = floor(dims ./ scale);  % nodes after scaling the original image in graph model
edges = prod(nodes) * neighbor; % edges for all nodes (each one has 8 or 4)
edgeCost = zeros(nodes(1), nodes(2), neighbor); % cost for each edge

tic;
intensity = conv2(img, ones(scale), 'valid');   % intensity matrix
intensity = intensity(1:scale:end, 1:scale:end);% each element is the sum of a group of pixels

for i = 1 : neighbor    % calculate the costs matrix for each neighbor
    edgeCost(:, :, i) = abs(intensity - imgshift(intensity, dn_area(i, :)));
end
%% generate MST model
nodeLabel = zeros(nodes);  % different tree set labels, 0 for single point out of trees
treeNodes = [];
treeEdges = [];
edgeCandidate = [];

for i = 1 : nodes(1)	% construct an edge list contains all edges in the node matrix
    for j = 1 : nodes(2)
        pos = [i, j];   % find a point
        pos_n = repmat(pos, [neighbor/2 1]) + fliplr(dn_area(1:neighbor/2, :)); % find it's 4 neighbors
        pos_n(:, end+1) = 1:neighbor/2;   % record edge's order
        [out1, ~] = find(pos_n(:, 1:end-1) < 1);    % find positions out of boundary
        [out2, ~] = find(pos_n(:, 1:end-1) > repmat(nodes, [neighbor/2 1]));
        outIndex = union(out1, out2);
        pos_n(outIndex, :) = [];    % remove positions outside intensity matrix
        list = [repmat(pos, [size(pos_n,1) 1]), pos_n(:, 1:end-1), ... % construct a small edge list
        reshape(edgeCost(pos(1), pos(2), pos_n(:, end)), [size(pos_n,1) 1])];
        edgeCandidate = [edgeCandidate; list]; % merge list with previous one
    end
end

[~, index] = sort(edgeCandidate(:, end)); % sort the edge candidates list in ascending order
edgeCandidate = edgeCandidate(index, :); % the smallest edge at the begining

set = 1;
% while(min(nodeLabel(:)) == 0 || numel(unique(nodeLabel)) > 1)
while(min(nodeLabel(:)) == 0)
    minEdge = []; 
    newNode = [];
    node1 = edgeCandidate(1, 1:2); node2 = edgeCandidate(1, 3:4); % get two nodes of the minimal edge
    label1 = nodeLabel(node1(1), node1(2)); label2 = nodeLabel(node2(1), node2(2)); % get labels of two nodes
    if label1 == 0 && label2 == 0   % two nodes are both not in trees
        minEdge = edgeCandidate(1, :);	% read minimal edge information
        newNode = [edgeCandidate(1, 1:2); edgeCandidate(1, 3:4)]; % read new node
        nodeLabel(node1(1), node1(2)) = set; nodeLabel(node2(1), node2(2)) = set; % assign new label to them
        set = set + 1; % label should be distinct for different tree sets
    elseif label1 * label2 == 0    % only one of two nodes is in tree
        minEdge = edgeCandidate(1, :);	% read minimal edge information
        if label1 == 0
            newNode = edgeCandidate(1, 1:2); % read new node
            nodeLabel(node1(1), node1(2)) = label2; % assign the label to one not labelled
        elseif label2 == 0
            newNode = edgeCandidate(1, 3:4); % read new node
            nodeLabel(node2(1), node2(2)) = label1; % assign the label to one not labelled
        end
    elseif label1 ~= label2    % two nodes are in two trees
        minEdge = edgeCandidate(1, :);      % read minimal edge information
        nodeLabel(nodeLabel == max(label1, label2)) = min(label1, label2); % union labels to the smaller one, two trees merge into one   
    elseif label1 == label2
        % two nodes are already in one tree, do nothing
    end
    
    edgeCandidate(1, :) = [];   % remove the minimal edge
    
    if ~isempty(newNode)    % find a new node with minimal edge
        treeNodes = [treeNodes; newNode];   % add the node into tree
    end
    if ~isempty(minEdge)   % find a valid edge
        treeEdges = [treeEdges; minEdge];   % add edge into edge set
    end
end

fprintf('MST calculation completes. ');
toc;
%% illustration for MST model
subplot(2, 2, 2);
imagesc(intensity); title(['downsampled image, ', ...
    num2str(nodes(1)), ' x ', num2str(nodes(2))]); axis off;

subplot(2, 2, 3);
drawMST(treeNodes, treeEdges, nodes, nodeLabel); title('generated MST'); 

[~, pos] = max(intensity(:)); 
recoveredImg = intensity .*  double(nodeLabel == nodeLabel(pos));
subplot(2, 2, 4);
imagesc(recoveredImg); title(['recovered image, ', num2str(nodes(1)), ...
    ' x ', num2str(nodes(2))]); axis off;
