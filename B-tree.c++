/*
Step1: Start Program w a tree notation, binary tree, how do we want to represent
Structure or array(keep track of node)

Step 2: Do any rotations on any edges

Goal: Given two sets of trees what’s the minimum of flips to get to object and how they are corresponding 

When a edge is common between beginning and end you don’t flip
 - Add a leaf for each node to make a complete tree, label leaf, traverse each node
 -can directly check if they are equivalent 

If we can determine a edge is same base on start and end without leafs

Write Question: Consider the binary rotation that allows me to quickly flip a edge and find similarities of edges between two trees

Brute Force as baseline, try to rotate them test from there
*/