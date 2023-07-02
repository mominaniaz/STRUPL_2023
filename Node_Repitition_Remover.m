
function nodal_connectivity_Repitition_Remover = Node_Repitition_Remover(nodal_connectivity_values)

[rows,columns] = size(nodal_connectivity_values);

nodal_connectivity_Repitition_Remover = zeros(rows,columns);
index_length = rows*columns;



nodal_connectivity_Repitition_Remover = nodal_connectivity_Repitition_Remover';
nodal_connectivity_values = nodal_connectivity_values';

nodal_connectivity_Repitition_Remover(1:3,1)  = nodal_connectivity_values(1:3,1);

check_if_same = false;

for index = 4:index_length
    check_if_same = false;
    for index_inside = 1:index-1
        if check_if_same
            continue;
        end
        if nodal_connectivity_values(index) == nodal_connectivity_values(index_inside)
            check_if_same = true;
            continue;
        end
    end
    if check_if_same == false
        nodal_connectivity_Repitition_Remover(index) = nodal_connectivity_values(index);
    end
    
end

nodal_connectivity_Repitition_Remover = nodal_connectivity_Repitition_Remover';

end

