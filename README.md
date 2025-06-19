# Current Objectives: 

Implement Brute Force Algorithm

**Progress so far**  
-- Implemented Helper Functions: TreeEquals, TreetoString, BFS(W/O Free edge logic)
-- Basic testing to ensure code is working as expected

**Todo List**  

-- Implement free edge shortcut
--Dist function


**Current Phase Question**  
   Can 2 edges of the same range in different trees be left and right edges? 
      - No, requires more evidence


**Brute Force Algorithm Approach**  

We compute FindRotationDistance(T₀,T₁) by calling

Dist(Tₛ,Tₑ):

If they’re identical, return 0.

If we’ve seen this pair (cached by their leaf‐range signatures), return the memo

If there’s a “free” edge e in Tₛ (rotating e inserts an edge already in Tₑ), do that rotation, split both trees along the new edge into two subtree pairs, and return the sum of Dist on each pair

Otherwise fall back to BFSSearch, which does a level‐by‐level BFS of all single rotations (left‐rotate/right‐rotate at each node), stopping as soon as you either match Tₑ or hit a free edge (and recurse immediately)

By memoizing on leaf‐range vectors and always peeling off free edges first (divide‐and‐conquer), with BFS only when necessary, you get an O*(cᵏ) FPT algorithm for rotation distance

**Paper Research**





