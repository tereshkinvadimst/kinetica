import medcoupling as mc


mesh: mc.MEDCouplingUMesh = mc.ReadUMeshFromFile("Mesh_1.med", "Mesh_1")


edgeMesh, desc, descI, revDesc, revDescI = mesh.buildDescendingConnectivity()

mesh.writeVTK("mesh.vtu")