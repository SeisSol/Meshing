#ifndef WRITER_H_
#define WRITER_H_

#include <cassert>
#include <fstream>
#include <iostream>

#include "hdf5.h"

#include "mesh.h"

enum class Hdf5DataType { hdf5UInt64, hdf5Int32, hdf5Double };

template <typename TT>
static void checkH5ErrImpl(TT status, const char* file, int line) {
  if (status < 0)
    std::cerr << "An HDF5 error occurred (" << file << ": " << line << ")" << std::endl;
}

#define checkH5Err(...) checkH5ErrImpl(__VA_ARGS__, __FILE__, __LINE__)

class Writer {
  public:
  explicit Writer(const std::string& fileNamePrefix, const Mesh& mesh)
      : mesh(mesh), xdmfFileName(fileNamePrefix + ".xdmf"), hdfFileName(fileNamePrefix + ".h5"){};

  void write(const std::vector<std::string>& parameterNames, const std::vector<std::vector<double>>& materialValues) {
    assert(parameterNames.size() == materialValues.size());
    writeXdmf(parameterNames);
    writeHdf5(parameterNames, materialValues);
  };

  private:
  const Mesh& mesh;
  const std::string xdmfFileName;
  const std::string hdfFileName;

  void writeHdf5(const std::vector<std::string>& parameterNames,
                 const std::vector<std::vector<double>>& materialValues);
  void writeXdmf(const std::vector<std::string>& parameterNames);
  void
      writeData(hid_t h5file, const std::string& name, std::vector<hsize_t> sizes, Hdf5DataType type, const void* data);
};
#endif
