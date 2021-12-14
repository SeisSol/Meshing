#ifndef PLACE_RECEIVERS_LINEARRECEIVERSEARCH_H
#define PLACE_RECEIVERS_LINEARRECEIVERSEARCH_H

#include <vector>

union Point {
  double coords[3];
  struct {
    double x;
    double y;
    double z;
  };
};

struct Triangle {
  Point vertices[3];
  double normals[3][2];
  double dist[3];
  double faceNormal[3];
  double faceDist;

  void determineNormals() {
    for (unsigned side = 0; side < 3; ++side) {
      double* v0 = vertices[side].coords;
      double* v1 = vertices[(side+1)%3].coords;
      double* vtest = vertices[(side+2)%3].coords;
      // normal = (y,-x)
      normals[side][0] = v1[1]-v0[1];
      normals[side][1] = -v1[0]+v0[0];

      dist[side] = normals[side][0] * v0[0] + normals[side][1] * v0[1];

      if (normals[side][0] * vtest[0] + normals[side][1] * vtest[1] > dist[side]) {
        normals[side][0] *= -1.;
        normals[side][1] *= -1.;
        dist[side] *= -1;
      }
    }

    double a[3];
    double b[3];
    for (unsigned d = 0; d < 3; ++d) {
      a[d] = vertices[1].coords[d]-vertices[0].coords[d];
      b[d] = vertices[2].coords[d]-vertices[0].coords[d];
    }
    faceNormal[0] = a[1]*b[2] - a[2]*b[1];
    faceNormal[1] = a[2]*b[0] - a[0]*b[2];
    faceNormal[2] = a[0]*b[1] - a[1]*b[0];
    faceDist = faceNormal[0] * vertices[0].coords[0] + faceNormal[1] * vertices[0].coords[1] + faceNormal[2] * vertices[0].coords[2];
  }

public:
  bool isInside(const Point& point) const {
    bool inside = true;
    for (unsigned side = 0; side < 3; ++side) {
      inside = inside && (normals[side][0] * point.x + normals[side][1] * point.y <= dist[side]);
    }
    return inside;
  }

  Triangle(Point vertex0, Point vertex1, Point vertex2) {
    vertices[0] = vertex0;
    vertices[1] = vertex1;
    vertices[2] = vertex2;
    determineNormals();
  }
};

class LinearReceiverSearch {
  std::vector<Point> receivers;
  std::vector<bool> processed;

public:
  LinearReceiverSearch(std::vector<Point> points) :
  receivers(std::move(points)) {
    processed = std::vector<bool>(receivers.size(), false);
  };

  template<typename Function>
  void applyFunctionOnAllReceiversInTriangle(const Triangle& triangle,
                                             Function f) {
      for (int i = 0; i < receivers.size(); ++i) {
        if (processed[i]) continue;
        auto& receiver = receivers[i];
        if (triangle.isInside(receiver)) {
          receiver = f(receiver);
          processed[i] = true;
        }
      }
  }

  const std::vector<Point>& getReceivers() const {
    return receivers;
  }
};

#endif // PLACE_RECEIVERS_LINEARRECEIVERSEARCH_H
