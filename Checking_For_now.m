
function nodal_connectivity_Repitition_Remover = Checking_For_now(nodal_connectivity_values);

[rows,columns] = size(nodal_connectivity_values);

nodal_connectivity_Repitition_Remover = zeros(rows,columns);


for i = 1:rows
    for j = 1:columns
        for i_t = 1:i
            for j_t = 1:j
                if nodal_connectivity_values(i,j) ~= nodal_connectivity_values(i_t,j_t)
                    nodal_connectivity_Repitition_Remover(i,j) = nodal_connectivity_values(i,j);
                end
            end
        end
    end
end


end

