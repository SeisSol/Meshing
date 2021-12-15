#ifndef READER_H_
#define READER_H_

#include <string>
#include <vector>
#include <optional>
#include <array>

class Mesh {
public:
    using neighbor_t = std::optional<unsigned>;
    using neighbors_t = std::array<neighbor_t, 4>;
    using boundaries_t = std::array<int, 4>;

    explicit Mesh(const std::string &fileName);

    unsigned int elementSize;
    unsigned int vertexSize;
    std::vector<boundaries_t> elementBoundaries;
    std::vector<int> elementGroups;
    std::vector<unsigned long> connect;
    std::vector<double> vertices;
    std::vector<neighbors_t> neighbors;

};

int decodeBoundaryCondition(int encoded, int faceId);

#endif
