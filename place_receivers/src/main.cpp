/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#include <cmath>
#include <iostream>
#include <utils/args.h>

#include "Reader.h"
#include "Geometry.h"

int main(int argc, char** argv)
{
  utils::Args args;
  args.addOption("receivers", 'r', "Receiver locations");
  args.addOption("mesh", 'm', "Netcdf mesh file");
  args.addOption("output", 'o', "Receiver output file");
	//~ args.addOption("mcs", 'm', "Proj.4 string that describes the mesh coordinate system (e.g. \"+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs\").", utils::Args::Required, false);
  //~ args.addOption("output", 'o', "Output file (.nrf)", utils::Args::Required, false);
  //~ args.addOption("normalize-onset", 'n', "Subtract the minimum onset time from all onsets.", utils::Args::No, false);
	//~ args.addOption("vcs", 'v', "Proj.4 string that describes the coordinate system for visualisation (defaults to geocentric, i.e. \"+proj=geocent +datum=WGS84 +units=m +no_def\").", utils::Args::Required, false);
  //~ args.addOption("xdmf", 'x', "Output for visualisation (.xmf)", utils::Args::Required, false);
  //~ 
  //~ args.setCustomHelpMessage("\nWith rconv you may either convert a SRF file to a NRF file, which you can use as input in SeisSol.\n"
                            //~ "In this case, give the options -i, -m, -o, and optionally -n.\n\n"
                            //~ "You may also write a file which may be loaded in Paraview for visualisation of the SRF file.\n"
                            //~ "In this case, give the options -i, -x, and optionally -v.\n\n"
                            //~ "You may write both files simultaneously by giving all options.\n");

	if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }
  
  std::string receiverFile = args.getArgument<std::string>("receivers");
  std::string meshFile = args.getArgument<std::string>("mesh");
  std::string receiverOutputFile = args.getArgument<std::string>("output");

  std::vector<Point> receivers = readReceiverFile(receiverFile);
  std::cout << "Read " << receivers.size() << " receivers." << std::endl;
  KDTree tree(receivers, 1);
  
  Mesh mesh(meshFile);
  
  for (unsigned p = 0; p < mesh.partitions; ++p) {
    if (p % 100 == 0) {
      std::cout << "Processing partition " << p << " to " << std::min(mesh.partitions, static_cast<size_t>(p+99)) << std::endl;
    }

    mesh.readPartition(p);
    
    setElevation(p, mesh, tree);
  }

  writeReceiverFile(tree, receiverOutputFile); 
  
  return 0;
}
