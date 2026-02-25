#ifndef REMOVE_WATERLAYER_WRITER_H
#define REMOVE_WATERLAYER_WRITER_H

#include <iostream>
#include <fstream>
#include "Reader.h"
#include "hdf5.h"

enum class Hdf5DataType {
    hdf5UInt64,
    hdf5Int32,
    hdf5Double
};

template<typename TT>
static void checkH5ErrImpl(TT status, const char *file, int line) {
  if (status < 0)
    std::cerr << "An HDF5 error occurred (" << file << ": " << line << ")" << std::endl;
}

#define checkH5Err(...) checkH5ErrImpl(__VA_ARGS__, __FILE__, __LINE__)

int encodeBoundary(std::array<int, 4> decoded) {
  int encodedBoundary = 0;
  for (auto i = 0u; i < 4; ++i) {
    encodedBoundary += decoded[i] << (i * 8u);
  }
  return encodedBoundary;
}

class Writer {
public:
    explicit Writer(Mesh *mesh) : mesh(mesh) {};

    void writeHdf5(const std::string &filename) {
      // Create file
      hid_t h5falist = H5Pcreate(H5P_FILE_ACCESS);
      checkH5Err(h5falist);
#ifdef H5F_LIBVER_V18
      checkH5Err(H5Pset_libver_bounds(h5falist, H5F_LIBVER_V18, H5F_LIBVER_V18));
#else
      checkH5Err(H5Pset_libver_bounds(h5falist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
#endif
      checkH5Err(H5Pset_fapl_mpio(h5falist, MPI_COMM_WORLD, MPI_INFO_NULL));
      hid_t h5file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5falist);
      checkH5Err(h5file);
      checkH5Err(H5Pclose(h5falist));

      writeData(h5file, "/connect",
                std::vector<hsize_t>{mesh->elementSize, 4},
                Hdf5DataType::hdf5UInt64,
                mesh->connect.data());

      writeData(h5file, "/geometry",
                std::vector<hsize_t>{mesh->vertexSize, 3},
                Hdf5DataType::hdf5Double,
                mesh->vertices.data());

      writeData(h5file, "/group",
                std::vector<hsize_t>({mesh->elementSize}),
                Hdf5DataType::hdf5Int32,
                mesh->elementGroups.data());

      auto encodedBoundary = std::vector<int>(mesh->elementSize);
      for (auto i = 0u; i < mesh->elementBoundaries.size(); ++i) {
        encodedBoundary[i] = encodeBoundary(mesh->elementBoundaries[i]);
      }
      writeData(h5file, "/boundary",
                std::vector<hsize_t>({mesh->elementSize}),
                Hdf5DataType::hdf5Int32,
                encodedBoundary.data());
    }

    void writeXdmf(const std::string &xdmfFileName, const std::string &h5FileName) {
      std::ofstream xdmf(xdmfFileName.c_str());

      const auto globalSize = std::array<unsigned int, 2>{mesh->elementSize, mesh->vertexSize};
      xdmf << "<?xml version=\"1.0\" ?>" << std::endl
           << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
           << "<Xdmf Version=\"2.0\">" << std::endl
           << " <Domain>" << std::endl
           << "  <Grid Name=\"puml mesh\" GridType=\"Uniform\">" << std::endl
           << "   <Topology TopologyType=\"Tetrahedron\" NumberOfElements=\"" << globalSize[0]
           << "\">"
           << std::endl
           // This should be UInt but for some reason this does not work with
           // binary data
           << "    <DataItem NumberType=\"Int\" Precision=\"8\" Format=\"HDF\" "
              "Dimensions=\""
           << globalSize[0] << " 4\">" << h5FileName << ":/connect</DataItem>" << std::endl
           << "   </Topology>" << std::endl
           << "   <Geometry name=\"geo\" GeometryType=\"XYZ\" NumberOfElements=\"" << globalSize[1]
           << "\">" << std::endl
           << "    <DataItem NumberType=\"Float\" Precision=\"" << sizeof(double)
           << "\" Format=\"HDF\" Dimensions=\"" << globalSize[1] << " 3\">" << h5FileName
           << ":/geometry</DataItem>" << std::endl
           << "   </Geometry>" << std::endl
           << "   <Attribute Name=\"group\" Center=\"Cell\">" << std::endl
           << "    <DataItem  NumberType=\"Int\" Precision=\"4\" Format=\"HDF\" "
              "Dimensions=\""
           << globalSize[0] << "\">" << h5FileName << ":/group</DataItem>" << std::endl
           << "   </Attribute>" << std::endl
           << "   <Attribute Name=\"boundary\" Center=\"Cell\">" << std::endl
           << "    <DataItem NumberType=\"Int\" Precision=\"4\" Format=\"HDF\" "
              "Dimensions=\""
           << globalSize[0] << "\">" << h5FileName << ":/boundary</DataItem>" << std::endl
           << "   </Attribute>" << std::endl
           << "  </Grid>" << std::endl
           << " </Domain>" << std::endl
           << "</Xdmf>" << std::endl;
    }

private:
    Mesh *mesh;

    void writeData(hid_t h5file,
                   const std::string &name,
                   std::vector<hsize_t> sizes,
                   Hdf5DataType type,
                   void *data) {
      hid_t h5space = H5Screate_simple(sizes.size(), sizes.data(), nullptr);
      checkH5Err(h5space);

      hid_t h5data;
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

};


#endif //REMOVE_WATERLAYER_WRITER_H
