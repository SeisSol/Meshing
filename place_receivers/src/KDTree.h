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

#ifndef TOOLS_KDTREE_H_
#define TOOLS_KDTREE_H_

#include <vector>

union Point {
  double coords[3];
  struct {
    double x;
    double y;
    double z;
  };
};

class KDTree {
public:
	KDTree(std::vector<Point> const& points, int maxLeafSize);
	~KDTree();
	
	template<typename Support, typename Action>
	void search(Support const& support, Action& action) const
	{
		searchTree<Support, Action>(0, support, action);
	}
	
	inline Point const* points() const { return data; }
	inline int index(int r) const { return idx[r]; }
  inline int numPoints() const { return nodes[0].n;}

private:
	int leftChild(int k) const { return 2*k + 1; }
	int rightChild(int k) const { return 2*k + 2; }

	struct Node {
    Node() : isLeaf(false) {}
		double pivot;
		int start;
		int n;
		int splitdim;
		bool isLeaf;
	};
	Node* nodes;
	
	void swap(int i, int j);
	int partition(int left, int right, int pivotIdx, int splitdim);
	void buildTree(int k, int splitdim);
	
	template<typename Support, typename Action>
	void searchTree(int k, Support const& support, Action& action) const;

	Point* data;
	int* idx;
	int p;
	int maxLeafN;
};

template<typename Support, typename Action>
void KDTree::searchTree(int k, Support const& support, Action& action) const
{
	Node& node = nodes[k];
	if (node.isLeaf) {
		for (int i = node.start; i < node.start + node.n; ++i) {
      action(data[i]);
		}
	} else {
		if (support(node.splitdim, 0) <= node.pivot) {
			searchTree<Support, Action>(leftChild(k), support, action);
		}
		if (support(node.splitdim, 1) >= node.pivot) {
			searchTree<Support, Action>(rightChild(k), support, action);			
		}
	}
}

#endif
