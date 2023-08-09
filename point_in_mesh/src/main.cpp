#include <iostream>
#include <utils/args.h>
#include <utils/logger.h>

#include "Geometry.h"
#include "Reader.h"

int main(int argc, char** argv) {
  utils::Args args;
  args.addOption("receivers", 'r', "SeisSol Receiver File");
  args.addOption("mesh", 'm', "Mesh File");

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  auto receiverFile = args.getArgument<std::string>("receivers");
  auto meshFile = args.getArgument<std::string>("mesh");

  auto receivers = reader::readReceiverFile(receiverFile);
  logInfo() << "Read" << receivers.size() << "receivers.";

  reader::Mesh const mesh(meshFile);

  geometry::PointChecker pointChecker(mesh);

  for (const auto& point : receivers) {
    const auto result = pointChecker.pointInMesh(point);
    if (result.has_value()) {
      logInfo() << "Found point (" << point[0] << ", " << point[1] << ", " << point[2]
                << ") in cell" << result.value() << ".";
    } else {
      logInfo() << "Did not find point (" << point[0] << ", " << point[1] << ", " << point[2]
                << ") "
                << "in any mesh cell.";
    }
  }

  return 0;
}
