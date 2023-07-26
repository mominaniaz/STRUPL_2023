function node_repitition = Node_Repitition_Counter(nodal_connectivity_values,number_of_nodes)

[rows,columns] = size(nodal_connectivity_values);
index_length = rows*columns;

%% Finding a matrix that will provide the number of times a node is repeated

node_repitition = zeros(number_of_nodes,1);

for i = 1:number_of_nodes
    times_repeated = 0;
    for index = 1:index_length
        if i == nodal_connectivity_values(index)
            times_repeated = times_repeated + 1;
        end
    end
    node_repitition(i) = times_repeated;
end

end