#include <cmath>
#include <iostream>
#include <utils/args.h>
#include <set>
#include <regex>
#include <iterator>
#include <algorithm>
#include <charconv>
#include <filesystem>

#include "mpi.h"
#include "Reader.h"
#include "Filter.h"
#include "Writer.h"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  utils::Args args;
  args.addOption("input", 'i', "input mesh");
  args.addOption("output", 'o', "output mesh");
  args.addOption("remove-groups", 'g', "groups to remove");

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }


  auto groupsString = args.getArgument<std::string>("remove-groups");
  const std::regex groupsRegex("\\,");
  std::set<unsigned int> groupsToRemove{};
  std::for_each(std::sregex_token_iterator(groupsString.begin(),
                                           groupsString.end(),
                                           groupsRegex,
                                           -1),
                std::sregex_token_iterator(),
                [&groupsToRemove](auto &arg) {
                    unsigned int groupInt;
                    auto[ptr, parseError] = std::from_chars(
                            arg.str().data(),
                            arg.str().data() + arg.str().size(),
                            groupInt);
                    if (parseError == std::errc::invalid_argument
                        || parseError == std::errc::result_out_of_range) {
                      std::cerr << "Group string is not a comma separated list of unsigned integers" << std::endl;
                      std::abort();
                    }
                    groupsToRemove.insert(groupInt);
                });


  // Input
  const auto meshFile = args.getArgument<std::string>("input");
  auto mesh = Mesh(meshFile);

  // Filtering
  removeGroups(mesh, groupsToRemove);

  // Output
  const auto outputFile = args.getArgument<std::filesystem::path>("output");
  auto xdmfFile = outputFile;
  xdmfFile.replace_extension("xdmf");

  Writer writer{&mesh};
  std::cout << "Writing hdf5 file: " << outputFile << std::endl;
  writer.writeHdf5(outputFile);
  std::cout << "Writing xdmf file: " << xdmfFile << std::endl;
  writer.writeXdmf(xdmfFile, "test.h5");

  return 0;
}
