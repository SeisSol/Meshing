/* (c) 2023 SeisSol, David Schneller. BSD-3 License */

#include <H5Fpublic.h>
#include <H5Spublic.h>
#include <mpi.h>
#include <hdf5.h>
#include <omp.h>

#include <array>
#include <string>
#include <vector>
#include <cmath>
#include <cstddef>

#include "utils/args.h"
#include "utils/logger.h"

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

        periodic = args.getArgument<bool>("bper" + dim2str[dimension], true);
        if (periodic) {
            // periodic boundary values
            minBoundary = 6;
            maxBoundary = 6;
        }
        else {
            minBoundary = args.getArgument<int>("bmin" + dim2str[dimension], 1);
            maxBoundary = args.getArgument<int>("bmax" + dim2str[dimension], 1);
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

    std::array<int, 6> cubeBoundary(std::size_t id) const {
        auto cid = unpackCubeId(id);
        std::array<int, 6> indices;
        for (int i = 0; i < 6; ++i) {
            indices[i] = 0;
            for (int j = 0; j < 4; ++j) {
                auto offset = cube2tetface[i][j];
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
                indices[i] |= ((type & 0xff) << (8*j));
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

    hid_t h5p = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_layout(h5p, H5D_CONTIGUOUS);
	H5Pset_alloc_time(h5p, H5D_ALLOC_TIME_EARLY);

    hid_t globalSpace = H5Screate_simple(globalDims.size(), globalDims.data(), nullptr);
    hid_t localSpace = H5Screate_simple(localDims.size(), localDims.data(), nullptr);
    hid_t datasetHandle = H5Dcreate(fileHandle, name.c_str(), datatype, globalSpace, H5P_DEFAULT, h5p, H5P_DEFAULT);

    H5Pclose(h5p);

    H5Sselect_hyperslab(globalSpace, H5S_SELECT_SET, localStart.data(), nullptr, localDims.data(), nullptr);
    H5Sselect_all(localSpace);

    std::vector<T> buffer(localFlattened);

    #pragma omp parallel for
    for (std::size_t i = 0; i < localCount; ++i) {
        std::invoke(accessor, i + localOffset, i * localChunk, buffer);
    }

    H5Dwrite(datasetHandle, datatype, localSpace, globalSpace, xfer, buffer.data());

    H5Sclose(localSpace);
    H5Sclose(globalSpace);

    H5Dclose(datasetHandle);
}

int main(int argc, char** argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int commSize;
    int commRank;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    utils::Args args;
    args.addOption("bperx", 0, "periodicity in x dimension", utils::Args::Required, false);
	args.addOption("bpery", 0, "periodicity in y dimension", utils::Args::Required, false);
	args.addOption("bperz", 0, "periodicity in z dimension", utils::Args::Required, false);
	args.addOption("bminx", 0, "boundary condition at the start of x dimension (ignored, if bperx is true)", utils::Args::Required, false);
	args.addOption("bmaxx", 0, "boundary condition at the end of x dimension (ignored, if bperx is true)", utils::Args::Required, false);
	args.addOption("bminy", 0, "boundary condition at the start of y dimension (ignored, if bpery is true)", utils::Args::Required, false);
	args.addOption("bmaxy", 0, "boundary condition at the end of y dimension (ignored, if bpery is true)", utils::Args::Required, false);
	args.addOption("bminz", 0, "boundary condition at the start of z dimension (ignored, if bperz is true)", utils::Args::Required, false);
	args.addOption("bmaxz", 0, "boundary condition at the end of z dimension (ignored, if bperz is true)", utils::Args::Required, false);
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

	logInfo(commRank) << "Using" << omp_get_max_threads() << "threads, and " << commSize << "ranks";

    if (args.parse(argc, argv)) {
        logError() << "Error parsing arguments";
    }

    Cubes cubes(args);

    auto filename = args.getArgument<std::string>("output");

    logInfo(commRank) << "Output to file" << filename;

    hid_t h5plist = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_libver_bounds(h5plist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
	H5Pset_fapl_mpio(h5plist, MPI_COMM_WORLD, MPI_INFO_NULL);

    hid_t fileHandle = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5plist);
    H5Pclose(h5plist);

    hid_t xfer = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);

    writeData<double>(fileHandle, xfer, "geometry", cubes.vertexCount(), 1, {3}, H5T_NATIVE_DOUBLE, [&](std::size_t id, std::size_t datapos, std::vector<double>& target) {
        auto data = cubes.vertex(id);
        target[datapos+0] = data[0];
        target[datapos+1] = data[1];
        target[datapos+2] = data[2];
    });
    writeData<unsigned long>(fileHandle, xfer, "connect", cubes.cubeCount(), 6, {4}, H5T_NATIVE_ULONG, [&](std::size_t id, std::size_t datapos, std::vector<unsigned long>& target){
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
    writeData<unsigned int>(fileHandle, xfer, "boundary", cubes.cubeCount(), 6, {}, H5T_NATIVE_UINT, [&](std::size_t id, std::size_t datapos, std::vector<unsigned int>& target){
        auto data = cubes.cubeBoundary(id);
        for (int i = 0; i < 6; ++i) {
            target[datapos + i] = data[i];
        }
    });
    writeData<unsigned long>(fileHandle, xfer, "identify", cubes.vertexCount(), 1, {}, H5T_NATIVE_ULONG, [&](std::size_t id, std::size_t datapos, std::vector<unsigned long>& target){
        target[datapos] = cubes.vertexIdentification(id);
    });

    H5Pclose(xfer);
    H5Fclose(fileHandle);

    MPI_Finalize();
    return 0;
}
