#include "writer.h"

void Writer::writeData(
    hid_t h5file, const std::string& name, const std::vector<hsize_t>& sizes, Hdf5DataType type, const void* data) {
  hid_t h5space = H5Screate_simple(sizes.size(), sizes.data(), nullptr);
  checkH5Err(h5space);

  hid_t h5data = 0;
  switch (type) {
  case Hdf5DataType::hdf5UInt64:
    h5data = H5Dcreate(h5file, name.c_str(), H5T_STD_U64LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    break;
  case Hdf5DataType::hdf5Int32:
    h5data = H5Dcreate(h5file, name.c_str(), H5T_STD_I32LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    break;
  case Hdf5DataType::hdf5Double:
    h5data = H5Dcreate(h5file, name.c_str(), H5T_IEEE_F64LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    break;
  }
  checkH5Err(h5data);

  hsize_t start[2] = {0, 0};
  checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, nullptr, sizes.data(), nullptr));

  hid_t h5memspace = H5Screate_simple(sizes.size(), sizes.data(), nullptr);
  checkH5Err(h5memspace);

  hid_t h5dxlist = H5Pcreate(H5P_DATASET_XFER);
  checkH5Err(h5dxlist);
  checkH5Err(H5Pset_dxpl_mpio(h5dxlist, H5FD_MPIO_COLLECTIVE));

  switch (type) {
  case Hdf5DataType::hdf5UInt64:
    checkH5Err(H5Dwrite(h5data, H5T_NATIVE_ULONG, h5memspace, h5space, h5dxlist, data));
    break;
  case Hdf5DataType::hdf5Int32:
    checkH5Err(H5Dwrite(h5data, H5T_NATIVE_INT, h5memspace, h5space, h5dxlist, data));
    break;
  case Hdf5DataType::hdf5Double:
    checkH5Err(H5Dwrite(h5data, H5T_NATIVE_DOUBLE, h5memspace, h5space, h5dxlist, data));
    break;
  }
}

void Writer::writeXdmf(const std::vector<std::string>& parameterNames) {
  std::ofstream xdmf(xdmfFileName.c_str());

  const auto nElements = mesh.getElementVertices().size();
  const auto nVertices = mesh.getVertexCoordinates().size();

  xdmf << R"(<?xml version="1.0" ?>)" << std::endl
       << R"(<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>)" << std::endl
       << R"(<Xdmf Version="2.0">)" << std::endl
       << R"( <Domain>)" << std::endl
       << R"(  <Grid Name="puml mesh" GridType="Uniform">)" << std::endl
       << R"(   <Topology TopologyType="Tetrahedron" NumberOfElements=")" << nElements << R"(">)"
       << std::endl
       // This should be UInt but for some reason this does not work with
       // binary data
       << R"(    <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions=")" << nElements << R"( 4">)"
       << hdfFileName << R"(:/connect</DataItem>)" << std::endl
       << R"(   </Topology>)" << std::endl
       << R"(   <Geometry name="geo" GeometryType="XYZ" NumberOfElements=")" << nVertices << R"(">)" << std::endl
       << R"(    <DataItem NumberType="Float" Precision=")" << sizeof(double) << R"(" Format="HDF" Dimensions=")"
       << nVertices << R"( 3">)" << hdfFileName << R"(:/geometry</DataItem>)" << std::endl
       << R"(   </Geometry>)" << std::endl
       << R"(   <Attribute Name="group" Center="Cell">)" << std::endl
       << R"(    <DataItem NumberType="Int" Precision="4" Format="HDF" Dimensions=")" << nElements << R"(">)"
       << hdfFileName << R"(:/group</DataItem>)" << std::endl
       << R"(   </Attribute>)" << std::endl;
  for (const auto& param : parameterNames) {
    xdmf << R"(   <Attribute Name=")" << param << R"(" Center="Cell">)" << std::endl
         << R"(    <DataItem  NumberType="Float" Precision=")" << sizeof(double) << R"(" Format="HDF" Dimensions=")"
         << nElements << R"(">)" << hdfFileName << R"(:/)" << param << R"(</DataItem>)" << std::endl
         << R"(   </Attribute>)" << std::endl;
  }
  xdmf << "  </Grid>" << std::endl << " </Domain>" << std::endl << "</Xdmf>" << std::endl;
}

void Writer::writeHdf5(const std::vector<std::string>& parameterNames,
                       const std::vector<std::vector<double>>& materialValues) {
  // Create file
  hid_t h5falist = H5Pcreate(H5P_FILE_ACCESS);
  checkH5Err(h5falist);
#ifdef H5F_LIBVER_V18
  checkH5Err(H5Pset_libver_bounds(h5falist, H5F_LIBVER_V18, H5F_LIBVER_V18));
#else
  checkH5Err(H5Pset_libver_bounds(h5falist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
#endif
  checkH5Err(H5Pset_fapl_mpio(h5falist, MPI_COMM_WORLD, MPI_INFO_NULL));
  hid_t h5file = H5Fcreate(hdfFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5falist);
  checkH5Err(h5file);
  checkH5Err(H5Pclose(h5falist));

  writeData(h5file,
            "/connect",
            std::vector<hsize_t>{mesh.getElementVertices().size(), 4},
            Hdf5DataType::hdf5UInt64,
            mesh.getElementVertices().data());

  writeData(h5file,
            "/geometry",
            std::vector<hsize_t>{mesh.getVertexCoordinates().size(), 3},
            Hdf5DataType::hdf5Double,
            mesh.getVertexCoordinates().data());

  writeData(h5file,
            "/group",
            std::vector<hsize_t>{mesh.getElementGroups().size()},
            Hdf5DataType::hdf5Int32,
            mesh.getElementGroups().data());

  for (unsigned i = 0; i < parameterNames.size(); i++) {
    writeData(h5file,
              "/" + parameterNames[i],
              std::vector<hsize_t>{mesh.getElementVertices().size()},
              Hdf5DataType::hdf5Double,
              materialValues[i].data());
  }
}
