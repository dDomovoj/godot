#ifndef OCTREE_TABLES_H
#define OCTREE_TABLES_H

namespace OctreeTables {

//  Corners:                                    Octants:
//
//         6---------------18--------------7       o---o---o
//        /               /               /|       | 6 | 7 |
//       /               /               / |       o---o---o  Upper
//      17--------------25--------------19 |       | 5 | 4 |
//     /               /               /   |       o---o---o
//    /               /               /    |
//   5---------------16--------------4     |       o---o---o
//   |     14--------|-----23--------|-----15      | 2 | 3 |
//   |    /          |    /          |    /|       o---o---o  Lower        Z
//   |   /           |   /           |   / |       | 1 | 0 |               |
//   |  22-----------|--26-----------|--24 |       o---o---o           X---o
//   | /             | /             | /   |
//   |/              |/              |/    |
//   13--------------21--------------12    |
//   |     2---------|-----10--------|-----3
//   |    /          |    /          |    /
//   |   /           |   /           |   /
//   |  9------------|--20-----------|--11           Y
//   | /             | /             | /             | Z
//   |/              |/              |/              |/
//   1---------------8---------------0          X----o

// The order is important for some algorithms (DMC)
const int g_octant_position[8][3]{

	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 1, 0, 1 },
	{ 0, 0, 1 },

	{ 0, 1, 0 },
	{ 1, 1, 0 },
	{ 1, 1, 1 },
	{ 0, 1, 1 }
};

} // namespace OctreeTables

#endif // OCTREE_TABLES_H
