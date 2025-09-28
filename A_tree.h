#ifndef A_TREE_H
#define A_TREE_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <utility>

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
    friend int  BFSSearch(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te);
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

std::unordered_set<RP,RPH> buildTargetSet(const VectorRangeTreeMap&);
bool hasEdgeByRange(const VectorRangeTreeMap& tree, const RP& e);
bool hasFreeEdge(const VectorRangeTreeMap& cur, const VectorRangeTreeMap& tgt,
                 int& out_v, bool& out_leftRotation);

// Distance Computation & Search
int BFSSearch(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te);
int removeFreeEdge(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te);
int Dist(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te);
int FindRotationDistance(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final);

// Verification tests
void runAllTests();

#endif // A_TREE_H
