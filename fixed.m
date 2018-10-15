function [R,K] = fixed(R,K,edge_array,T)
edge_nodes = [edge_array(:,1)' edge_array(end,end)];
    for i = 1:length(edge_nodes)
        R(edge_nodes(i)) = T;
        K(edge_nodes(i),:) = sparse(1,length(K));
        K(edge_nodes(i),edge_nodes(i)) = 1;
    end
    
end