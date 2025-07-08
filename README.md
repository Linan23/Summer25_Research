# Current Objectives: 

Implement Brute Force Algorithm

**Progress so far**  
-- Implemented Helper Functions: TreeEquals, TreetoString, BFS(W/O Free edge logic)
-- Basic testing to ensure code is working as expected
-- Organizing file into cpp & header for better readability, accessability, optimization
-- Implemented all basic functions of provided sudo code, with unit testing to verify correctness

**Todo List**  

Some other detail functions to keep in mind....

-- Scaling/Profiling 
-- Optimize Range updates: perhaps incrementally update only affected path
-- Harden Free-Edge Detection
-- Memoization optimization??
-- Additional Testing(eg empty trees,single nodes, skewed etc )



**Current Phase Question**  
   Can 2 edges of the same range in different trees be left and right edges? 
      - No


**Brute Force Algorithm Approach**  

We compute FindRotationDistance(T₀,T₁) by calling

Dist(Tₛ,Tₑ):

If they’re identical, return 0.

If we’ve seen this pair (cached by their leaf‐range signatures), return the memo

If there’s a “free” edge e in Tₛ (rotating e inserts an edge already in Tₑ), do that rotation, split both trees along the new edge into two subtree pairs, and return the sum of Dist on each pair

Otherwise fall back to BFSSearch, which does a level‐by‐level BFS of all single rotations (left‐rotate/right‐rotate at each node), stopping as soon as you either match Tₑ or hit a free edge (and recurse immediately)

By memoizing on leaf‐range vectors and always peeling off free edges first (divide‐and‐conquer), with BFS only when necessary, you get an O*(cᵏ) FPT algorithm for rotation distance

**Paper Algorithm**





