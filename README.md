## Underactuated Consensus Between Manipulator and Mobile Robot Thesis
This repository includes the MATLAB simulation code of my thesis project. The abstract of the thesis is as follows:

### Abstract:
Robotics has been a huge field of study over the past few decades. Since generally robots are mechanical structures capable of moving, aiming to assist humans in their tasks, we can meet various types of robots, such as ground mobile robots, aerial robots (Unmanned Aerial Vehicle - UAV), underwater robots, robotic arms (manipulators), or even combinations of those types, such as mobile manipulators.

For our purposes, we will focus on the modeling and control of ground-vehicle mobile robots and robotic manipulators. Applications of mobile robotics can be found anywhere from simple warehouse tasks to self-driving cars and unknown outdoor environment exploration. On the other hand, robotic manipulators have mainly industrial usage.

This work is dedicated to the consensus between a ground mobile robot and a robotic manipulator so that they meet to a common location. More specifically, the goal is to control the two robots in a way that they reach a common pose with some given constraints while trying to maximize a given function concerning the manipulator. 

There are in general two approaches for the problem, the offline approach and the online approach, which are both investigated in this document. The offline approach consists of two distinct steps, the first step is to a-priori determine a suitable final position for the two robots to meet and the second step is to control each robot independently so that they both reach their previously selected pose. Meanwhile, the online approach involves a direct control of the robots, which attempts to achieve our goal, without any previous knowledge of the robots' final destination.

In this thesis, there are proposed some control schemes for both the offline and the online case, based mostly on the redundancy analysis regarding the manipulator and on the control algorithms for unicycles, regarding the mobile robot. The proposed solutions are firstly analyzed mathematically, in order to secure the control scheme's stability and afterward, they are applied in simulations so that their validity is confirmed.

---
---

Each of the files `Lagrange.m`, `Inclined_Redundant.m` and `Online.m` implement a different simulation for different approaches of the same problem. The full document for the project can be found at this [link](https://drive.google.com/file/d/1rDmU60s3x-fGvnMoR9PkPPZQP0YzS23R/view?usp=sharing), while a more detailed explanation of the code can be found in the appendix of the same document.

