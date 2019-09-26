% adapted from: 
% https://www.mathworks.com/matlabcentral/fileexchange/11947-b-a-scale-free-network-generation-and-visualization
% the original codes might generate self-loop
% Yu Zhu, Rice ECE, 03/28/2019

% -Nodes is the desired network size
% -mlinks controls the number of links a new node can make to the existing network nodes
% -seed is the original network to which the B-A algorithm links additional nodes with a specific preferential attachment procedure. 
% This undirected adjacency matrix can be created manually, or one could use the Adjacency Matrix GUI. 
% Each node must have at least one link. The seed variable can be replaced with a developed scale-free network to generate a larger one. 
% Make sure the new Nodes variable is greater than the size of the seed network.

function [A,L] = generate_connected_BA(Nodes, mlinks, seed)

  %  seed = full(seed);
    flag_connected = 0;

    while flag_connected ==0

        pos = length(seed);

        %if (Nodes < pos) || (mlinks > pos) || (pos < 1) || (size(size(seed)) ~= 2) || (mlinks < 1) || (seed ~= seed') || (sum(diag(seed)) ~= 0)
        %    error('invalid parameter value(s)');
        %end

        %if mlinks > 5 || Nodes > 15000 || pos > 15000
        %    warning('Abnormally large value(s) may cause long processing time');
        %end

        rand('state',sum(100*clock));

       % Net = zeros(Nodes, Nodes, 'single');
        Net = zeros(Nodes, Nodes);
        Net(1:pos,1:pos) = seed;
        sumlinks = sum(sum(Net));

        while pos < Nodes
            pos = pos + 1;
            linkage = 0;
            while linkage ~= mlinks
                rnode = ceil(rand * (pos-1)); % error fixed here
                deg = sum(Net(:,rnode)) * 2;
                rlink = rand * 1;
                if rlink < deg / sumlinks && Net(pos,rnode) ~= 1 && Net(rnode,pos) ~= 1
                    Net(pos,rnode) = 1;
                    Net(rnode,pos) = 1;
                    linkage = linkage + 1;
                    sumlinks = sumlinks + 2;
                end
            end
        end

        %clear Nodes deg linkage pos rlink rnode sumlinks mlinks
        clear deg linkage pos rlink rnode sumlinks 
        A = Net;
        L = diag(sum(A)) - A; 
        if sum(abs(eig(L)) <= 1.00e-06) == 1
            flag_connected = 1; % Check whether the graph is connected
        end
        
    end

end

% This file is used to simulate the B-A algorithm and returns scale-free networks of given node sizes. 
% Understanding the B-A algorithm is key to using this code to its fullest. 
% Due to Matlab resource limitations, it may not be possible to generate networks much larger than 15000 nodes, 
% and increasing the mlinks variable increases processing time severely. 
% This code was developed so that one could generate a network of small size, 
% and then use that network as a seed to build a greater sized network, 
% continuing this process until the actual desired network size is reached. 
% This is for processor and time management purposes. 
% However, realize that the initial seed does not have to have scale-free properties, 
% while the later seeds may happen to have these properties. 
% Therefore, it is prudent not to make the initial seed size much larger than a few nodes (most commonly 5 interconnected nodes). 
% In addition, the mlinks should be kept constant throughout the creation of the scale-free network.