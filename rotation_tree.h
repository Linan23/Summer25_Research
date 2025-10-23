// Defines the `VectorRangeTreeMap` structure—the shared binary tree
// representation used by both solvers and the comparison tooling. The class
// stores traversal metadata, supports in-place rotations, and exposes helpers
// for hashing, equality, and edge collection. Search code in bfs.cpp and the
// Java bridge depend on this header for consistent tree semantics.

#ifndef ROTATION_TREE_H
#define ROTATION_TREE_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <utility>

// Forward decls
struct BFSStats;

// Utility for hashing (parent, child) pairs
struct PairHash { size_t operator()(const std::pair<int,int>& p) const; };
struct PairEq   { bool operator()(const std::pair<int,int>& a, const std::pair<int,int>& b) const; };


// Core data structure
struct VectorRangeTreeMap {
    // data members
    std::vector<std::pair<int,int>> ranges;
    std::vector<std::pair<int,int>> edges;
    std::vector<int> parents;
    int root;
    int max_node_value;
    std::vector<int> original_inorder, original_preorder;
    std::unordered_map<int,int> position_in_inorder;
    std::unordered_set<int> original_nodes;

    static constexpr int NO_CHILD  = -1;
    static constexpr int NO_PARENT = -1;

    // ctor
    VectorRangeTreeMap();
    void recomputeUpwardsFrom(int start);

    // build
    void build(const std::vector<int>& preorder, const std::vector<int>& inorder);

#ifndef NDEBUG
    void verify() const;   // heavy checks (debug only)
#else
    inline void verify() const {}
#endif

    // accessors
    int getLeftChild(int node_value)  const;
    int getRightChild(int node_value) const;
    int getParent(int node_value)     const;
    std::pair<int,int> getRange(int node_value) const;
    bool isOriginal(int node_value) const;

    // mutators
    void setLeftChild(int node_value, int child_value);
    void setRightChild(int node_value, int child_value);

    void rotateLeft(int x);
    void rotateRight(int x);

    // partition
    static std::pair<VectorRangeTreeMap, VectorRangeTreeMap>
    partitionAlongEdge(const VectorRangeTreeMap& T,
                       const std::pair<int,int>& parent_range,
                       const std::pair<int,int>& child_range);

    void collectEdges(int node,
        std::unordered_set<std::pair<int,int>,PairHash,PairEq>& out) const;

    // friends
    friend bool        TreesEqual(const VectorRangeTreeMap& A, const VectorRangeTreeMap& B);
    friend std::string treeToString(const VectorRangeTreeMap& T);
    friend std::string canonicalTraversalString(const VectorRangeTreeMap& T);
    friend int  BFSSearch(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te,
                          BFSStats* stats);
    friend int  removeFreeEdge(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te);
    friend int  Dist(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te);
    friend int  FindRotationDistance(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final);

private:
    void buildRecursive(const std::vector<int>& preorder, const std::vector<int>& inorder,
                        int ps, int pe, int is, int ie, int parent);
    void clear();
    void calculateAllRanges();
    void updateNodeRange(int node_value);
    void calculateRangesPostOrder(int node_value);
};

// Free-Edge Detection Helpers
struct RP { int ps, pe, cs, ce; bool operator==(RP const&) const; };
struct RPH { size_t operator()(RP const&) const; };

static uint64_t fnv1a64(const std::string& s);

// Captures the outcome and telemetry of BFSSearchCapped.
struct BFSRun {
    int dist; size_t expanded, enqueued, visited;
    bool timeout, cap_hit; double seconds;
    size_t duplicates{0};
    size_t max_queue{0};
    size_t equality_hits{0};
    size_t signature_mismatches{0};
};

// Feature-rich BFS wrapper that enforces optional caps and returns telemetry.
BFSRun BFSSearchCapped(const VectorRangeTreeMap& T_s,
                       const VectorRangeTreeMap& T_e,
                       double time_limit_sec = 10.0,
                       size_t visited_cap    = 5'000'000,
                       size_t queue_cap      = 5'000'000);

// Optional live statistics callers can provide to BFSSearch.
struct BFSStats {
    size_t generated{0};
    size_t enqueued{0};
    size_t duplicates{0};
    size_t equality_hits{0};
    size_t max_queue{0};
    size_t visited{0};
    size_t signature_mismatches{0};
};

std::unordered_set<RP,RPH> buildTargetSet(const VectorRangeTreeMap&);
bool hasEdgeByRange(const VectorRangeTreeMap& tree, const RP& e);
bool hasFreeEdge(const VectorRangeTreeMap& cur, const VectorRangeTreeMap& tgt,
                 int& out_v, bool& out_leftRotation);

// Distance Computation & Search
// Full BFS with optional statistics (used by testing/CLI code).
int BFSSearch(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te,
              BFSStats* stats = nullptr);
// Baseline breadth-first search without heuristics.
int BFSSearchBaseline(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te,
                      BFSStats* stats = nullptr);
// Enhanced breadth-first search with hashing/filters/meet-in-the-middle helpers.
int BFSSearchOptimized(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te,
                       BFSStats* stats = nullptr);
// Hashed bidirectional BFS (meet-in-the-middle) with optional pruning heuristics.
int BiBFSSearchHashed(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te,
                      size_t state_cap = 2'000'000);
// Attempts to resolve a “free edge” via decomposition; INT_MAX means fallback.
int removeFreeEdge(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te);
// Memoised distance computation leveraging free-edge shortcuts plus BFS.
int Dist(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te);
// Public API used by clients to compute rotation distance end-to-end.
int FindRotationDistance(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final);

// Verification tests
void runAllTests();

// Serialisation helpers shared across CLI/test code.
std::string treeToString(const VectorRangeTreeMap& T);
std::string canonicalTraversalString(const VectorRangeTreeMap& T);

#endif // ROTATION_TREE_H
