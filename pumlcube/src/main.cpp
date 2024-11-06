// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include <hdf5.h>
#include <functional>
#include <mpi.h>
#include <hdf5.h>
#include <omp.h>

#include <array>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <stack>
#include <cassert>

#include "utils/args.h"
#include "utils/logger.h"

template <typename TT> static TT _hw(TT&& status, const char* file, int line) {
  if (status < 0) {
    logError() << utils::nospace << "An HDF5 error occurred (" << file << ": " << line << ")";
  }
  return std::forward<TT>(status);
}

#define hw(...) _hw(__VA_ARGS__, __FILE__, __LINE__)

enum class BoundaryFormat {
    Int32,
    Int64,
    Int32x4
};

// cf. https://en.wikipedia.org/wiki/Schl%C3%A4fli_orthoscheme
static std::array<std::array<std::array<int, 3>, 4>, 6> cube2tet = {
    std::array<std::array<int, 3>, 4>{
        std::array<int, 3>{0,0,0},
        std::array<int, 3>{1,0,0},
        std::array<int, 3>{0,1,0},
        std::array<int, 3>{0,1,1},
    },
    std::array<std::array<int, 3>, 4>{
        std::array<int, 3>{0,0,0},
        std::array<int, 3>{1,0,0},
        std::array<int, 3>{0,1,1},
        std::array<int, 3>{0,0,1},
    },
    std::array<std::array<int, 3>, 4>{
        std::array<int, 3>{1,0,0},
        std::array<int, 3>{0,1,1},
        std::array<int, 3>{0,0,1},
        std::array<int, 3>{1,0,1},
    },
    std::array<std::array<int, 3>, 4>{
        std::array<int, 3>{1,1,1},
        std::array<int, 3>{0,1,1},
        std::array<int, 3>{1,0,1},
        std::array<int, 3>{1,0,0},
    },
    std::array<std::array<int, 3>, 4>{
        std::array<int, 3>{1,1,1},
        std::array<int, 3>{0,1,1},
        std::array<int, 3>{1,0,0},
        std::array<int, 3>{1,1,0},
    },
    std::array<std::array<int, 3>, 4>{
        std::array<int, 3>{0,1,1},
        std::array<int, 3>{1,0,0},
        std::array<int, 3>{1,1,0},
        std::array<int, 3>{0,1,0},
    },
};

// I don't know why this order works. But it seems like to does.
static std::array<int, 4> faceorder = { 3,2,0,1 };

static std::array<std::array<std::array<int, 3>, 4>, 6> tet2face(const std::array<std::array<std::array<int, 3>, 4>, 6>& data) {
    std::array<std::array<std::array<int, 3>, 4>, 6> output;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int c = 0; c < 3; ++c) {
                std::array<int, 3> values;
                int d = 0;
                for (int k = 0; k < 4; ++k) {
                    if (faceorder[j] != k) {
                        values[d] = data[i][k][c];
                        ++d;
                    }
                }
                if (values[0] == values[1] && values[1] == values[2]) {
                    output[i][j][c] = values[0];
                }
                else {
                    output[i][j][c] = -1;
                }
            }
        }
    }
    return output;
}

static std::array<std::array<std::array<int, 3>, 4>, 6> cube2tetface = tet2face(cube2tet);

static std::array<std::string, 3> dim2str{"x", "y", "z"};

struct CubeDimension {
    bool periodic;
    int minBoundary, maxBoundary;
    std::size_t numCells;
    double scaling;
    double translation;

    // (NOTE: utils::Args does not support const& right now for what we want to do)
    CubeDimension(utils::Args& args, int dimension) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        logInfo(rank) << "Reading dimension" << dim2str[dimension];

        minBoundary = args.getArgument<int>("bmin" + dim2str[dimension], 1);
        maxBoundary = args.getArgument<int>("bmax" + dim2str[dimension], 1);

        periodic = minBoundary == 6 && maxBoundary == 6;
        if ((minBoundary == 6 || maxBoundary == 6) && !periodic) {
            logError() << "One-sided periodic boundary conditions (i.e. value 6 for bmin or bmax) are not supported.";
        }

        numCells = args.getArgument<std::size_t>("n" + dim2str[dimension]);

        scaling = args.getArgument<double>("s" + dim2str[dimension], 1.0);
        translation = args.getArgument<double>("t" + dim2str[dimension], 0.0);
    }

    std::size_t cellCount() const {
        return numCells;
    }

    std::size_t faceCount() const {
        return cellCount() + 1;
    }

    int getFaceType(std::size_t id) const {
        if (id == 0) {
            // left boundary
            return minBoundary;
        }
        else if (id == cellCount()) {
            // right boundary
            return maxBoundary;
        }
        else {
            // interior faces
            return 0;
        }
    }

    double facePosition(std::size_t id) const {
        return (scaling / cellCount()) * id + translation;
    }

    std::size_t getFaceIdentification(std::size_t id) const {
        // for periodicity, identify with the opposite-side vertex
        if (periodic && id == cellCount()) {
            return 0;
        }
        
        // otherwise, identify with itself
        return id;
    }
};

struct Cubes {
    CubeDimension dimx;
    CubeDimension dimy;
    CubeDimension dimz;

    Cubes(utils::Args& args)
        : dimx(args, 0), dimy(args, 1), dimz(args, 2)
    {
        
    }

    std::size_t vertexCount() const {
        return dimx.faceCount() * dimy.faceCount() * dimz.faceCount();
    }

    std::size_t cubeCount() const {
        return dimx.cellCount() * dimy.cellCount() * dimz.cellCount();
    }

    std::size_t tetrahedronCount() const {
        return cubeCount() * 6;
    }

    std::array<std::size_t, 3> unpackCubeId(std::size_t id) const {
        auto z = id / (dimy.cellCount() * dimx.cellCount());
        auto xy = id % (dimy.cellCount() * dimx.cellCount());
        auto y = xy / dimx.cellCount();
        auto x = xy % dimx.cellCount();
        return {x,y,z};
    }

    std::array<std::size_t, 3> unpackVertexId(std::size_t id) const {
        auto z = id / (dimy.faceCount() * dimx.faceCount());
        auto xy = id % (dimy.faceCount() * dimx.faceCount());
        auto y = xy / dimx.faceCount();
        auto x = xy % dimx.faceCount();
        return {x,y,z};
    }

    std::size_t packCubeId(const std::array<std::size_t, 3>& id) const {
        return dimx.cellCount() * dimy.cellCount() * id[2] + dimx.cellCount() * id[1] + id[0];
    }

    std::size_t packVertexId(const std::array<std::size_t, 3>& id) const {
        return dimx.faceCount() * dimy.faceCount() * id[2] + dimx.faceCount() * id[1] + id[0];
    }

    std::array<std::array<std::size_t, 4>, 6> cubeVertices(std::size_t id) const {
        auto cid = unpackCubeId(id);
        std::array<std::array<std::size_t, 4>, 6> indices;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 4; ++j) {
                auto offset = cube2tet[i][j];
                std::array<std::size_t, 3> vid = {cid[0] + offset[0], cid[1] + offset[1], cid[2] + offset[2]};
                indices[i][j] = packVertexId(vid);
            }
        }
        return indices;
    }

    int faceBoundary(std::size_t id, int i, int j) const {
        const auto cid = unpackCubeId(id);
        const auto offset = cube2tetface[i][j];
        int type = 0;
        if (offset[0] >= 0) {
            type = dimx.getFaceType(cid[0] + offset[0]);
        }
        if (offset[1] >= 0) {
            type = dimy.getFaceType(cid[1] + offset[1]);
        }
        if (offset[2] >= 0) {
            type = dimz.getFaceType(cid[2] + offset[2]);
        }
        return type;
    }

    template<typename T>
    std::array<T, 6> cubeBoundary(std::size_t id, int typesize) const {
        std::array<T, 6> indices;
        const T mask = (T(1) << typesize) - 1;
        for (int i = 0; i < 6; ++i) {
            indices[i] = 0;
            for (int j = 0; j < 4; ++j) {
                const auto type = faceBoundary(id, i, j);
                indices[i] |= ((type & mask) << (typesize*j));
            }
        }
        return indices;
    }

    std::array<int, 24> cubeBoundary4x(std::size_t id) const {
        std::array<int, 24> indices;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 4; ++j) {
                const auto type = faceBoundary(id, i, j);
                indices[i * 4 + j] = type;
            }
        }
        return indices;
    }

    std::array<double, 3> vertex(std::size_t id) const {
        auto vid = unpackVertexId(id);
        return {dimx.facePosition(vid[0]), dimy.facePosition(vid[1]), dimz.facePosition(vid[2])};
    }

    std::size_t vertexIdentification(std::size_t id) const {
        auto vid = unpackVertexId(id);
        return packVertexId({dimx.getFaceIdentification(vid[0]), dimy.getFaceIdentification(vid[1]), dimz.getFaceIdentification(vid[2])});
    }
};

class XmlWriter {
public:
    XmlWriter(const std::string& file) : stream(file) {

    }

    void writeIndent() {
        for (std::size_t i = 0; i < nodes.size() * indent; ++i) {
            stream << ' ';
        }
    }

    void writeHeader() {
        stream << "<?xml version=\"1.0\"?>\n";
    }

    void writeDoctype() {
        stream << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    }

    void writeNode(const std::string& name) {
        writeNodeHeaderEnd();
        writeIndent();
        stream << "<" << name;
        nodes.push(name);
        inHeader = true;
    }

    void writeNodeEnd() {
        const auto node = nodes.top();
        nodes.pop();
        if (inHeader) {
            stream << " />\n";
        }
        else {
            writeIndent();
            stream << "</" << node << ">\n";
        }
    }

    void writeNodeHeaderEnd() {
        if (inHeader) {
            stream << ">\n";
            inHeader = false;
        }
    }

    template<typename T>
    void writeNodeAttribute(const std::string& name, const T& value) {
        assert(inHeader);
        stream << " " << name << "=\"" << value << "\"";
    }

    template<typename T>
    void writeData(const T& data) {
        writeNodeHeaderEnd();
        writeIndent();
        stream << data << "\n";
    }

    void closeDocument() {
        while (!nodes.empty()) {
            writeNodeEnd();
        }
        stream.flush();
    }

    void writeComment(const std::string& text) {
        writeNodeHeaderEnd();
        writeIndent();
        stream << "<!--" << text << "-->\n";
    }
private:
    std::ofstream stream;
    std::stack<std::string> nodes;
    bool inHeader{false};
    int indent{2};
};

class XdmfXml {
public:
    XdmfXml(const std::string& fileName, const std::string& dataFile) : xml(fileName), dataFile(dataFile) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        localWrite = rank == 0;

        logInfo(rank) << "Writing Xdmf file to" << fileName;
    }

    void begin() {
        if (localWrite) {
            xml.writeHeader();
            xml.writeDoctype();
            xml.writeNode("Xdmf");
            xml.writeNodeAttribute("Version", "2.0");
            xml.writeNode("Domain");
            xml.writeNode("Grid");
            xml.writeNodeAttribute("Name", "cubemesh");
            xml.writeNodeAttribute("GridType", "Uniform");
        }
    }

    void addData(const std::string& dataset, const std::vector<std::size_t>& dimensions, const std::string& type, const std::string& name, const std::string& datatype, std::size_t size) {
        if (localWrite) {
            std::ostringstream dimensionStream;
            dimensionStream << dimensions[0];
            for (std::size_t i = 1; i < dimensions.size(); ++i) {
                dimensionStream << " " << dimensions[i];
            }
            if (type == "Topology") {
                xml.writeNode("Topology");
                xml.writeNodeAttribute("TopologyType", "Tetrahedron");
                xml.writeNodeAttribute("NumberOfElements", dimensions[0]);
            }
            else if (type == "Geometry") {
                xml.writeNode("Geometry");
                xml.writeNodeAttribute("GeometryType", "XYZ");
                xml.writeNodeAttribute("NumberOfElements", dimensions[0]);
            }
            else if (type == "AttributeCell") {
                xml.writeNode("Attribute");
                xml.writeNodeAttribute("Name", name);
                xml.writeNodeAttribute("Center", "Cell");
            }
            else if (type == "AttributeNode") {
                xml.writeNode("Attribute");
                xml.writeNodeAttribute("Name", name);
                xml.writeNodeAttribute("Center", "Node");
            }

            xml.writeNode("DataItem");
            xml.writeNodeAttribute("NumberType", datatype);
            xml.writeNodeAttribute("Precision", size);
            xml.writeNodeAttribute("Format", "HDF");
            xml.writeNodeAttribute("Dimensions", dimensionStream.str());
            xml.writeData(dataFile + ":/" + dataset);
            xml.writeNodeEnd();
            xml.writeNodeEnd();
        }
    }

    void end() {
        if (localWrite) {
            xml.closeDocument();
        }
    }
private:
    XmlWriter xml;
    bool localWrite;
    std::string dataFile;
};

template<typename T, typename F>
static void writeData(hid_t fileHandle, hid_t xfer, const std::string& name, std::size_t count, std::size_t chunksize, const std::vector<std::size_t>& elemsize, hid_t datatype, F&& accessor) {
    int commSize;
    int commRank;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    logInfo(commRank) << "Writing" << name;

    std::size_t localCount = count / commSize;
    std::size_t localOffset = localCount * commRank;
    auto rest = count % commSize;
    if ((std::size_t)commRank < rest) {
        ++localCount;
        localOffset += commRank;
    }
    else {
        localOffset += rest;
    }

    std::vector<hsize_t> globalDims(elemsize.size() + 1);
    globalDims[0] = chunksize * count;
    for (std::size_t i = 0; i < elemsize.size(); ++i) {
        globalDims[i + 1] = elemsize[i];
    }

    std::size_t localFlattened = 0;
    std::size_t localChunk = 0;
    std::vector<hsize_t> localDims(elemsize.size() + 1);
    localDims[0] = chunksize * localCount;
    localFlattened = localDims[0];
    localChunk = chunksize;
    for (std::size_t i = 0; i < elemsize.size(); ++i) {
        localDims[i + 1] = elemsize[i];
        localFlattened *= elemsize[i];
        localChunk *= elemsize[i];
    }

    std::vector<hsize_t> localStart(elemsize.size() + 1);
    localStart[0] = chunksize * localOffset;
    for (std::size_t i = 0; i < elemsize.size(); ++i) {
        localStart[i + 1] = 0;
    }

    hid_t h5p = hw(H5Pcreate(H5P_DATASET_CREATE));
    hw(H5Pset_layout(h5p, H5D_CONTIGUOUS));
    hw(H5Pset_alloc_time(h5p, H5D_ALLOC_TIME_EARLY));

    hid_t globalSpace = hw(H5Screate_simple(globalDims.size(), globalDims.data(), nullptr));
    hid_t localSpace = hw(H5Screate_simple(localDims.size(), localDims.data(), nullptr));
    hid_t datasetHandle = hw(H5Dcreate(fileHandle, name.c_str(), datatype, globalSpace, H5P_DEFAULT, h5p, H5P_DEFAULT));

    hw(H5Pclose(h5p));

    hw(H5Sselect_hyperslab(globalSpace, H5S_SELECT_SET, localStart.data(), nullptr, localDims.data(), nullptr));
    hw(H5Sselect_all(localSpace));

    std::vector<T> buffer(localFlattened);

    #pragma omp parallel for
    for (std::size_t i = 0; i < localCount; ++i) {
        std::invoke(accessor, i + localOffset, i * localChunk, buffer);
    }

    hw(H5Dwrite(datasetHandle, datatype, localSpace, globalSpace, xfer, buffer.data()));

    hw(H5Sclose(localSpace));
    hw(H5Sclose(globalSpace));

    hw(H5Dclose(datasetHandle));
}

int main(int argc, char** argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int commSize;
    int commRank;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    utils::Args args;
    args.addOption("bminx", 0, "boundary condition at the start of x dimension", utils::Args::Required, false);
    args.addOption("bmaxx", 0, "boundary condition at the end of x dimension", utils::Args::Required, false);
    args.addOption("bminy", 0, "boundary condition at the start of y dimension", utils::Args::Required, false);
    args.addOption("bmaxy", 0, "boundary condition at the end of y dimension", utils::Args::Required, false);
    args.addOption("bminz", 0, "boundary condition at the start of z dimension", utils::Args::Required, false);
    args.addOption("bmaxz", 0, "boundary condition at the end of z dimension", utils::Args::Required, false);
    args.addOption("nx", 'x', "number of cubes in x dimension", utils::Args::Required, true);
    args.addOption("ny", 'y', "number of cubes in y dimension", utils::Args::Required, true);
    args.addOption("nz", 'z', "number of cubes in z dimension", utils::Args::Required, true);
    args.addOption("output", 'o', "output file for resulting PUML mesh", utils::Args::Required, true);
    args.addOption("sx", 0, "scale in x direction", utils::Args::Required, false);
    args.addOption("sy", 0, "scale in y direction", utils::Args::Required, false);
    args.addOption("sz", 0, "scale in z direction", utils::Args::Required, false);
    args.addOption("tx", 0, "mesh translation in x direction (after scaling)", utils::Args::Required, false);
    args.addOption("ty", 0, "mesh translation in y direction (after scaling)", utils::Args::Required, false);
    args.addOption("tz", 0, "mesh translation in z direction (after scaling)", utils::Args::Required, false);
    args.addOption("xdmf", 0, "create an XDMF XML file", utils::Args::Required, false);
    args.addOption("boundary", 0, "boundary format to use", utils::Args::Required, false);

    logInfo(commRank) << "Using" << omp_get_max_threads() << "threads, and " << commSize << "ranks";

    if (args.parse(argc, argv)) {
        logError() << "Error parsing arguments";
    }

    Cubes cubes(args);

    const auto filename = args.getArgument<std::string>("output");

    const auto boundaryRaw = args.getArgument<std::string>("boundary", "i32");

    BoundaryFormat boundaryFormat = BoundaryFormat::Int32;
    if (boundaryRaw == "i32") {
        boundaryFormat = BoundaryFormat::Int32;
    }
    else if (boundaryRaw == "i64") {
        boundaryFormat = BoundaryFormat::Int64;
    }
    else if (boundaryRaw == "i32x4") {
        boundaryFormat = BoundaryFormat::Int32x4;
    }
    else {
        logError() << "Unknown boundary format:" << boundaryRaw;
    }

    const bool periodic = cubes.dimx.periodic || cubes.dimy.periodic || cubes.dimz.periodic;

    logInfo(commRank) << "Output to file" << filename;

    hid_t h5plist = hw(H5Pcreate(H5P_FILE_ACCESS));
    hw(H5Pset_libver_bounds(h5plist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
    hw(H5Pset_fapl_mpio(h5plist, MPI_COMM_WORLD, MPI_INFO_NULL));

    hid_t fileHandle = hw(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5plist));
    hw(H5Pclose(h5plist));

    hid_t xfer = hw(H5Pcreate(H5P_DATASET_XFER));
    hw(H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE));

    writeData<double>(fileHandle, xfer, "geometry", cubes.vertexCount(), 1, {3}, H5T_NATIVE_DOUBLE, [&](std::size_t id, std::size_t datapos, std::vector<double>& target) {
        auto data = cubes.vertex(id);
        target[datapos+0] = data[0];
        target[datapos+1] = data[1];
        target[datapos+2] = data[2];
    });
    writeData<uint64_t>(fileHandle, xfer, "connect", cubes.cubeCount(), 6, {4}, H5T_NATIVE_UINT64, [&](std::size_t id, std::size_t datapos, std::vector<uint64_t>& target){
        auto data = cubes.cubeVertices(id);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 4; ++j) {
                target[datapos + 4*i + j] = data[i][j];
            }
        }
    });
    writeData<int>(fileHandle, xfer, "group", cubes.cubeCount(), 6, {}, H5T_NATIVE_INT, [&](std::size_t id, std::size_t datapos, std::vector<int>& target){
        for (int i = 0; i < 6; ++i) {
            target[datapos + i] = 0;
        }
    });
    if (boundaryFormat == BoundaryFormat::Int32) {
        writeData<int32_t>(fileHandle, xfer, "boundary", cubes.cubeCount(), 6, {}, H5T_NATIVE_INT32, [&](std::size_t id, std::size_t datapos, std::vector<int32_t>& target){
            auto data = cubes.cubeBoundary<int>(id, 8);
            for (int i = 0; i < 6; ++i) {
                target[datapos + i] = data[i];
            }
        });
    }
    if (boundaryFormat == BoundaryFormat::Int64) {
        writeData<int64_t>(fileHandle, xfer, "boundary", cubes.cubeCount(), 6, {}, H5T_NATIVE_INT64, [&](std::size_t id, std::size_t datapos, std::vector<int64_t>& target){
            auto data = cubes.cubeBoundary<int64_t>(id, 16);
            for (int i = 0; i < 6; ++i) {
                target[datapos + i] = data[i];
            }
        });
    }
    if (boundaryFormat == BoundaryFormat::Int32x4) {
        writeData<int>(fileHandle, xfer, "boundary", cubes.cubeCount(), 6, {4}, H5T_NATIVE_INT, [&](std::size_t id, std::size_t datapos, std::vector<int>& target){
            auto data = cubes.cubeBoundary4x(id);
            for (int i = 0; i < 24; ++i) {
                target[datapos + i] = data[i];
            }
        });
    }
    if (periodic) {
        writeData<uint64_t>(fileHandle, xfer, "identify", cubes.vertexCount(), 1, {}, H5T_NATIVE_UINT64, [&](std::size_t id, std::size_t datapos, std::vector<uint64_t>& target){
            target[datapos] = cubes.vertexIdentification(id);
        });
    }

    hid_t attrSpace = hw(H5Screate(H5S_SCALAR));
    hid_t attrType = hw(H5Tcopy(H5T_C_S1));
    hw(H5Tset_size(attrType, H5T_VARIABLE));
    hid_t attrBoundary = hw(
        H5Acreate(fileHandle, "boundary-format", attrType, attrSpace, H5P_DEFAULT, H5P_DEFAULT));
    const void* stringData = boundaryRaw.data();
    hw(H5Awrite(attrBoundary, attrType, &stringData));
    hw(H5Aclose(attrBoundary));
    hw(H5Sclose(attrSpace));
    hw(H5Tclose(attrType));

    if (args.getArgument<bool>("xdmf", false)) {
        auto xdmf = XdmfXml(filename + ".xdmf", filename);
        xdmf.begin();
        xdmf.addData("geometry", {cubes.vertexCount(), 3}, "Geometry", "", "Float", 8);
        xdmf.addData("connect", {cubes.tetrahedronCount(), 4}, "Topology", "", "Int", 8);
        xdmf.addData("group", {cubes.tetrahedronCount()}, "AttributeCell", "group", "Int", 4);
        if (boundaryFormat == BoundaryFormat::Int32) {
            xdmf.addData("boundary", {cubes.tetrahedronCount()}, "AttributeCell", "boundary", "Int", 4);
        }
        if (boundaryFormat == BoundaryFormat::Int64) {
            xdmf.addData("boundary", {cubes.tetrahedronCount()}, "AttributeCell", "boundary", "Int", 8);
        }
        if (boundaryFormat == BoundaryFormat::Int32x4) {
            xdmf.addData("boundary", {cubes.tetrahedronCount(), 4}, "AttributeCell", "boundary", "Int", 4);
        }
        if (periodic) {
            xdmf.addData("identify", {cubes.vertexCount()}, "AttributeNode", "identify", "Int", 8);
        }
        xdmf.end();
    }

    hw(H5Pclose(xfer));
    hw(H5Fclose(fileHandle));

    logInfo(commRank) << "All done.";

    MPI_Finalize();
    return 0;
}
