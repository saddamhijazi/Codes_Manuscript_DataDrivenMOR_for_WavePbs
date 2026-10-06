function B = constructB(nmb_nodes, nmb_eNodes)
% Construct the matrix B with which one can add the Dirichlet boundary conditions
    B = sparse(nmb_nodes, nmb_nodes-nmb_eNodes);
    for j = 1:nmb_nodes-nmb_eNodes
        B(nmb_eNodes+j,j) = 1;
    end
end
