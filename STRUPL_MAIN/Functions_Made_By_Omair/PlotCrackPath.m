function PlotCrackPath(Crack_Path,nodal_coordinte_values,nodal_connectivity_values,Direction)

%% Drawings the Graph and the Elements with Nodes
nel = length(nodal_connectivity_values);
nnd = length(nodal_coordinte_values);
nne = size(nodal_connectivity_values,2);

X = zeros(nne,nel);
Y = zeros(nne,nel);

for iel=1:nel
    for i=1:nne
        nd(i)=nodal_connectivity_values(iel,i);         % extract connected node for (iel)-th element
        X(i,iel)=nodal_coordinte_values(nd(i),1);    % extract x value of the node
        Y(i,iel)=nodal_coordinte_values(nd(i),2);    % extract y value of the node
    end
end

    % Plotting the FEM mesh, diaplay Node numbers and Element numbers
    f1 = figure ;
    set(f1,'name','Mesh','numbertitle','off') ;
    plot(X,Y,'k')
    fill(X,Y,'w')
    
    title('Finite Element Mesh') ;
    axis off ;

%% Utilizing the Crack_path_matrix to draw the crack paths

Crack_Path_Matrix_Max = sortrows(Crack_Path,2);

node_configuration = Node_Configurator(nodal_coordinte_values);

node_configuration = node_configuration';

[rows,columns] = size(Crack_Path_Matrix_Max);

%% Using the highest Cracks from Crack_Path_matrix_max and making paths to the 
% nearest horizontal vertical etc 
for i = 1:rows

end







end