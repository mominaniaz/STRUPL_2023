function PlotCrackPath(nodal_coordinte_values,nodal_connectivity_values)

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

end