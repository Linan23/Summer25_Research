# Current Objectives: 

**Tree Representation**  
   - Decide how to store a binary tree:  
     - As a linked structure (nodes with left/right pointers)  
     - Or as an array (with parent/child indices)  

**Edge Rotations**  
   - Apply rotations to any edge in the tree.  

**Overall Goal**  
   - Given two binary trees, find the minimum number of flips (rotations) to turn one into the other.  
   - Track which nodes/edges correspond between the start and end trees.  

**Skip Common Edges**  
   - If an edge exists in both trees, don’t rotate it.  
   - Two ways to detect common edges:  
     1. **Complete the trees**  
        - Add dummy leaves so every internal node has two children.  
        - Label leaves uniquely, traverse, and compare edges directly.  
     2. **Direct check**  
        - Devise a method to spot matching edges without adding leaves.  

**Current Phase Question**  
   > How can we use binary rotations to quickly flip edges and identify shared edges between two trees?  

**Baseline Approach(ATM)**  
   - Start with a brute-force rotation tester
   - Measure its performance and use that as a reference before moving to FPT optimizations