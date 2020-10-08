/*
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Thomas Ulrich 
 *
 * @section LICENSE
 * Copyright (c) 2014-2015, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 
19.02.2016
Create a planar fault of given dimension and dip angle, the nucleation is also included in the geometry
To use obtain the mesh:
gmsh planar-anydip.geo -3 -optimize
To convert the mesh:
gmsh2gambit -i planar-anydip.msh -o planar-anydip.neu

 */


mesh_size = 10e3;
Xmax = -160e3;
Xmin = 215e3;
Ymin = 1235e3;
Ymax = 1605e3;
Zmin=-200e3;
Zmax=5e3;

//Create the Volume
Point(1) = {Xmin, Ymin, Zmax, mesh_size};
Point(2) = {Xmin, Ymax, Zmax, mesh_size};
Point(3) = {Xmax, Ymax, Zmax, mesh_size};
Point(4) = {Xmax, Ymin, Zmax, mesh_size};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1,2,3,4};
Plane Surface(1) = {5};

Extrude {0,0, Zmin} { Surface{1}; }

Physical Surface(101) = {1,14,18,22,26,27};
//Physical Surface(103) = {100,200};
//This ones are read in the gui
//Physical Surface(105) = {14,18,22,26,27};

//Physical Volume(1) = {1};
