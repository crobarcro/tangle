// SPDX-License-Identifier: MIT
// Copyright (c) 2026 David Meeker

// tangle_mesh.h — Tangle mesh structures and library API
#ifndef TANGLE_MESH_H
#define TANGLE_MESH_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <array>
#include <utility>

// Minimal tangle struct definitions needed by the bridge code.
// These must match the definitions in tangle.cpp.

struct Point {
    double x, y;
    int id;
    int marker;
    std::vector<double> attribs;
    double lfs = -1.0;
};

struct Triangle {
    std::array<int,3> v;
    std::array<int,3> neighbors;
    double region_attrib;
    double region_max_area;
    int marker;
    int generation = 0;
};

struct Segment {
    int v0, v1;
    int marker;
    double lfs = -1.0;
    int pbc_type = -1;
    bool no_split = false;
    int pbc_side = -1;   // 0/1: declared side of an (anti)periodic pair; -1 = n/a
};

struct Hole { double x, y; };
struct Region { double x, y; double attrib; double max_area; };

struct PBCNodePair {
    int node_a, node_b;
    int type;
};

struct PBCDef {
    // A PBC declaration from the input: the boundary with marker_a pairs with
    // the boundary with marker_b (.poly form; FEMM inputs share one marker).
    int marker_a, marker_b;
    int type; // 0=periodic, 1=anti-periodic
};

struct AGEDef {
    std::string name;
    int format;
    double innerAngle, outerAngle;
    double ri, ro;
    double totalArcLength;
    double cx, cy;
    int n;
    std::vector<int> innerNodes;
    std::vector<int> outerNodes;
};

// Ordered correspondence between two matched boundary chains, independent of
// any periodic field semantics. nodes_a[i] and nodes_b[i] are the final mesh
// nodes that refinement kept in one-to-one correspondence; marker_a/marker_b
// record which source boundary each chain came from. This is emitted whenever
// a PBC declaration pairs two chains, whether or not the caller wants the
// periodic node pairs that make up pbc_pairs.
struct BoundaryChainMatch {
    int marker_a;
    int marker_b;
    int type; // 0=periodic, 1=anti-periodic
    std::vector<int> nodes_a;
    std::vector<int> nodes_b;
};

struct Mesh {
    std::vector<Point>    vertices;
    std::vector<Triangle> triangles;
    std::vector<Segment>  segments;
    std::vector<Hole>     holes;
    std::vector<Region>   regions;
    std::vector<std::pair<int,int>> edges;
    std::vector<char> edge_boundary;     // parallel to edges: 1 if boundary (1 incident triangle)
    std::vector<PBCNodePair> pbc_pairs;
    std::map<int,int> pbc_twin;
    std::map<int,int> pbc_node_type;
    std::vector<PBCDef> pbc_defs;
    std::vector<AGEDef> age_defs;
    std::vector<BoundaryChainMatch> boundary_matches; // ordered matched chains

    void rebuildAdjacency();
    int locateTriangle(double px, double py, int hint = 0) const;
    bool hasEdge(int a, int b, const std::vector<std::vector<int>>& v2t) const;
};

// Library/process return status. Each distinct failure egress has its own code
// so a caller (e.g. femm-qt's "Mesher failed (error N)") can tell a missing file
// from a parse error from a real meshing failure. Keep in sync with the copy in
// tangle.cpp.
enum TangleStatus {
    TANGLE_OK          = 0,  // success
    TANGLE_ERR_USAGE   = 1,  // bad invocation: no args / no input file named
    TANGLE_ERR_NO_FILE = 2,  // input file not found / could not be opened
    TANGLE_ERR_PARSE   = 3,  // file opened but could not be parsed
    TANGLE_ERR_MESH    = 4,  // meshing failed (e.g. a segment could not be enforced)
    TANGLE_ERR_OPTION  = 5,  // option not valid for this input (e.g. -g on a non-FEMM file)
};

// Optional library-side overrides for FEMM meshing. A zero/false value leaves
// the value derived from the input problem unchanged.
struct MeshOptions {
    double minimumAngleDegrees = 0.0;       // 0 = use the problem's [MinAngle]
    double maximumElementArea = 0.0;        // >0 = apply to every region
    bool   forceMaximumElementArea = false; // true = override larger region limits
    bool   suppressExteriorSteinerPoints = false;
    bool   suppressUnusedVertices = false;
    bool   verbose = false;
};

// In-memory FEMM-like problem description accepted by the record-based library
// entry point. Field conventions match the .fem file columns: boundary and
// circuit references are 1-based, and 0 means "none". This lets a caller mesh
// a problem it has already parsed without a file round trip, while sharing the
// exact file reader (arc discretization, LFS, PBC and AGE handling).
struct FemProblem {
    bool   isMagnetics = true;  // .fem magnetics numbering (4/5/6/7 BdryFormat)
    double minAngle = 20.0;
    bool   doSmartMesh = false;

    struct Node {
        double x = 0.0, y = 0.0;
        int boundaryMarker = 0;  // 1-based; 0 = none
        int group = 0;
    };
    struct Segment {
        int n0 = 0, n1 = 0;
        double maxSideLength = 0.0;
        int boundaryMarker = 0;  // 1-based; 0 = none
        int hidden = 0;
        int group = 0;
    };
    struct Arc {
        int n0 = 0, n1 = 0;
        double arcLength = 0.0;
        double maxSegDegrees = 0.0;
        int boundaryMarker = 0;  // 1-based; 0 = none
        int hidden = 0;
        int group = 0;
    };
    struct Label {
        double x = 0.0, y = 0.0;
        int blockType = 0;             // 1-based; 0 = <None>
        double maxAreaDiameter = -1.0; // .fem mesh-size column (diameter); <=0 = none
        int inCircuit = 0;             // 1-based; 0 = none
        double magDir = 0.0;
        int group = 0;
        int turns = 1;
        int isExternal = 0;
    };
    struct Boundary {
        std::string name;
        int format = 0;                // 4/5 periodic/antiperiodic, 6/7 AGE
        double innerAngle = 0.0;
        double outerAngle = 0.0;
    };

    std::vector<Node> nodes;
    std::vector<Segment> segments;
    std::vector<Arc> arcs;
    std::vector<Label> labels;
    std::vector<Boundary> boundaries;
};

// Mesh a .fem file in-memory. Returns TANGLE_OK (0) on success, else one of the
// TANGLE_ERR_* codes above (tangle_mesh_fem yields NO_FILE / PARSE / MESH).
// The MeshOptions overload applies the caller overrides before meshing; the
// no-options overload uses the problem-derived defaults.
int tangle_mesh_fem(const std::string& inputBase, Mesh& outMesh);
int tangle_mesh_fem(const std::string& inputBase, const MeshOptions& options,
                    Mesh& outMesh);

// Mesh an in-memory FEMM-like problem. Same semantics and status codes as the
// file overload, but reads no file. Additive: the file overloads are unchanged.
int tangle_mesh_fem(const FemProblem& problem, const MeshOptions& options,
                    Mesh& outMesh);

#endif
