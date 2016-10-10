#include "KDTree.h"

#include <algorithm>
#include <cstring>
#include <cmath>

KDTree::KDTree(std::vector<Point> const& points, int maxLeafSize)
	: p(2), maxLeafN(maxLeafSize)
{	
	int n = points.size();
	
	// Copy data and change row-major to column-major storage.
	data = new Point[n];
	memcpy(data, points.data(), n * sizeof(Point));
	
	idx = new int[n];
	for (int i = 0; i < n; ++i) {
		idx[i] = i;
	}

	int maxHeight = 1 + ceil(log2(n / static_cast<double>(maxLeafN)));
	int maxNodes = (1 << maxHeight) - 1;
	nodes = new Node[maxNodes];
	nodes[0].start = 0;
	nodes[0].n = n;
	
	buildTree(0, 0);
}

KDTree::~KDTree()
{
	delete[] data;
	delete[] idx;
	delete[] nodes;
}

void KDTree::swap(int i, int j)
{
	if (i != j) {
		std::swap(idx[i], idx[j]);
    std::swap(data[i], data[j]);
	}
}

int KDTree::partition(int left, int right, int pivotIdx, int splitdim)
{
	double pivot = data[pivotIdx].coords[splitdim];
	int st = left;
	for (int i = left; i < right; ++i) {
		if (data[i].coords[splitdim] < pivot) {
			swap(st, i);
			++st;
		}
	}
	swap(st, right);
	return st;
}

void KDTree::buildTree(int k, int splitdim)
{
	Node& node = nodes[k];
	node.splitdim = splitdim;
	if (node.n > maxLeafN) {
		int half = (node.n % 2 != 0) ? (node.n + 1)/2 : node.n/2;
		int l = node.start;
		int r = l + node.n - 1;
		int median_idx = l + half;
		if (l != r) {
			int pivotIdx;
			while (true) {
				pivotIdx = r;
				pivotIdx = partition(l, r, pivotIdx, splitdim);
				if (median_idx == pivotIdx) {
					break;
				} else if (median_idx < pivotIdx) {
					r = pivotIdx - 1;
				} else {
					l = pivotIdx + 1;
				}
			}
		}
		node.pivot = data[median_idx].coords[splitdim];

		Node& left = nodes[leftChild(k)];
		left.start = node.start;
		left.n = half;

		Node& right = nodes[rightChild(k)];
		right.start = median_idx;
		right.n = node.n - half;
		
		int nextSplitdim = (splitdim + 1) % p;
		buildTree(leftChild(k), nextSplitdim);
		buildTree(rightChild(k), nextSplitdim);
	} else {
		node.isLeaf = true;
	}
}
