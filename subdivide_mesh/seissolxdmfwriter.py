import h5py
import numpy as np


def write_mesh_xdmf(
    prefix, nNodes, nCells, lDataName, lData, node_per_element, reduce_precision
):
    precisionDict = {"int64": 8, "int32": 4, "float64": 8, "float32": 4}
    numberTypeDict = {
        "int64": "UInt",
        "int32": "UInt",
        "float64": "Float",
        "float32": "Float",
    }
    topology = "Tetrahedron" if node_per_element == 4 else "Triangle"
    xdmf = f"""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
  <Grid Name="puml mesh" GridType="Uniform">
   <Topology TopologyType="{topology}" NumberOfElements="{nCells}">
    <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="{nCells} {node_per_element}">{prefix}.h5:/connect</DataItem>
   </Topology>
   <Geometry name="geo" GeometryType="XYZ" NumberOfElements="{nNodes}">
    <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="{nNodes} 3">{prefix}.h5:/geometry</DataItem>
   </Geometry>"""
    for k, dataName in enumerate(lDataName):
        prec = 4 if reduce_precision else precisionDict[lData[k].dtype.name]
        number_type = numberTypeDict[lData[k].dtype.name]
        xdmf += f"""
    <Attribute Name="{dataName}" Center="Cell">
      <DataItem NumberType="{number_type}" Precision="{prec}" Format="HDF" Dimensions="1 {nCells}">{prefix}.h5:/{dataName}</DataItem>
    </Attribute>"""
    xdmf += """
  </Grid>
 </Domain>
</Xdmf>
"""
    with open(prefix + ".xdmf", "w") as fid:
        fid.write(xdmf)
    print(f"done writing {prefix}.xdmf")


def write_timeseries_xdmf(
    prefix,
    nNodes,
    nCells,
    lDataName,
    lData,
    dt,
    node_per_element,
    lidt,
    reduce_precision,
):
    precisionDict = {"int64": 8, "int32": 4, "float64": 8, "float32": 4}
    numberTypeDict = {
        "int64": "UInt",
        "int32": "UInt",
        "float64": "Float",
        "float32": "Float",
    }
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
        for k, dataName in enumerate(lDataName):
            prec = 4 if reduce_precision else precisionDict[lData[k].dtype.name]
            number_type = numberTypeDict[lData[k].dtype.name]
            xdmf += f"""
    <Attribute Name="{dataName}" Center="Cell">
     <DataItem ItemType="HyperSlab" Dimensions="{nCells}">
      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">{i} 0 1 1 1 {nCells}</DataItem>
      <DataItem NumberType="{number_type}" Precision="{prec}" Format="HDF" Dimensions="{i+1} {nCells}">{prefix}.h5:/{dataName}</DataItem>
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


def write_h5(prefix, lDataName, xyz, connect, lData, lidt, reduce_precision):
    dtypeDict = {
        "int64": "i8",
        "int32": "i4",
        "float64": "float64",
        "float32": "float32",
    }
    reducePrecisionDict = {
        "i8": "i4",
        "i4": "i4",
        "float64": "float32",
        "float32": "float32",
    }
    nCells, node_per_element = connect.shape
    with h5py.File(prefix + ".h5", "w") as h5f:
        h5f.create_dataset("/connect", (nCells, node_per_element), dtype="uint64")
        h5f["/connect"][:, :] = connect[:, :]
        h5f.create_dataset("/geometry", xyz.shape, dtype="d")
        h5f["/geometry"][:, :] = xyz[:, :]
        for k, dataName in enumerate(lDataName):
            hdname = "/" + dataName
            mydtype = dtypeDict[lData[k].dtype.name]
            if reduce_precision:
                mydtype = reducePrecisionDict[mydtype]
            if len(lData[0].shape) == 1 and len(lidt) == 1:
                h5f.create_dataset(hdname, (nCells,), dtype=str(mydtype))
                h5f[hdname][:] = lData[k][:]
            else:
                h5f.create_dataset(hdname, (len(lidt), nCells), dtype=str(mydtype))
                for i, idt in enumerate(lidt):
                    h5f[hdname][i, :] = lData[k][idt, :]
    print(f"done writing {prefix}.h5")


def write_seissol_output(
    prefix, xyz, connect, lDataName, lData, dt, lidt, reduce_precision=False
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
    reduce_precision: convert double to float and i64 to i32 if True
    """
    nNodes = xyz.shape[0]
    nCells, node_per_element = connect.shape
    if len(lidt) == 1:
        write_mesh_xdmf(
            prefix, nNodes, nCells, lDataName, lData, node_per_element, reduce_precision
        )
    else:
        write_timeseries_xdmf(
            prefix,
            nNodes,
            nCells,
            lDataName,
            lData,
            dt,
            node_per_element,
            lidt,
            reduce_precision,
        )
    write_h5(prefix, lDataName, xyz, connect, lData, lidt, reduce_precision)
