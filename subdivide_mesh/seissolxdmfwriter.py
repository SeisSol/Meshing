import h5py
import numpy as np

def write_seissol_xdmf(
    prefix, nNodes, nCells, lDataName, dt, node_per_element, lidt, precision
):
    prec = 8 if precision == "double" else 4
    topology = "Tetrahedron" if node_per_element == 4 else "Triangle"
    xdmf = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>"""
    for i, idt in enumerate(lidt):
        xdmf += f"""
  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
   <Grid Name="step_{idt}" GridType="Uniform">
    <Topology TopologyType="{topology}" NumberOfElements="{nCells}">
     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="{nCells} {node_per_element}">{prefix}.h5:/connect</DataItem>
    </Topology>
    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="{nNodes}">
     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="{nNodes} 3">{prefix}.h5:/geometry</DataItem>
    </Geometry>
    <Time Value="{idt*dt}"/>"""
        for dataName in lDataName:
            xdmf += f"""
    <Attribute Name="{dataName}" Center="Cell">
     <DataItem ItemType="HyperSlab" Dimensions="{nCells}">
      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">{i} 0 1 1 1 {nCells}</DataItem>
      <DataItem NumberType="Float" Precision="{prec}" Format="HDF" Dimensions="{i+1} {nCells}">{prefix}.h5:/{dataName}</DataItem>
     </DataItem>
    </Attribute>"""
        xdmf += """
   </Grid>
  </Grid>"""
    xdmf += """
 </Domain>
</Xdmf>
"""
    with open(prefix + ".xdmf", "w") as fid:
        fid.write(xdmf)
    print(f"done writing {prefix}.xdmf")


def write_seissol_h5(prefix, lDataName, xyz, connect, lData, lidt, precision):
    dtypeDict = {"int64": "i8", "int32": "i4", "float64": "float64", "float32": "float32"}
    nCells, node_per_element = connect.shape
    with h5py.File(prefix + ".h5", "w") as h5f:
        h5f.create_dataset("/connect", (nCells, node_per_element), dtype="uint64")
        h5f["/connect"][:, :] = connect[:, :]
        h5f.create_dataset("/geometry", xyz.shape, dtype="d")
        h5f["/geometry"][:, :] = xyz[:, :]
        for k, dataName in enumerate(lDataName):
            hdname = "/" + dataName
            h5f.create_dataset(hdname, (nCells), dtype=dtypeDict[lData[k].dtype.name])
            if len(lData[0].shape) == 1 and len(lidt) == 1:
                h5f[hdname][:] = lData[k][:]
            else:
                for i, idt in enumerate(lidt):
                    h5f[hdname][i, :] = lData[k][idt, :]
    print(f"done writing {prefix}.h5")

def write_seissol_output(
    prefix, xyz, connect, lDataName, lData, dt, lidt, precision="double"
):
    """
    Write hdf5/xdmf files output, readable by ParaView using SeisSol data
    prefix: file
    xyz: geometry array
    connect: connect array
    lDataName: list of array e.g. ['ASl', 'Vr']
    lData: list of numpy data array
    dt: sampling time of output
    lidt: list of time steps to be written
    """
    nNodes = xyz.shape[0]
    nCells, node_per_element = connect.shape
    write_seissol_xdmf(
        prefix, nNodes, nCells, lDataName, dt, node_per_element, lidt, precision
    )
    write_seissol_h5(prefix, lDataName, xyz, connect, lData, lidt, precision)
