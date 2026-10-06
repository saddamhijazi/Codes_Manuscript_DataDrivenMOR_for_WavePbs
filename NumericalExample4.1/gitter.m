function [nodes,mesh,model] = gitter(dx)

    numberOfPDE = 1;
    model = createpde(numberOfPDE);

    R1 = [3,4,0.0,5.0,5.0,0.0,5.0,5.0,0.0,0.0]';
    g = decsg(R1);
    geometryFromEdges(model,g);
    mesh = generateMesh(model,"GeometricOrder","linear", "Hmax",dx);
    nodes = model.Mesh.Nodes;

end 