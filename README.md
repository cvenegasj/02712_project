# Simulation of Cellular Dynamics with the Cellular Potts Model

Final project for the Course 02-712: Computational Methods for Biological Modeling and Simulation from Carnegie Mellon 
University, Fall 2019.

For this project, we aim to demonstrate the versatility of the Cell Potts Model (CPM) by extending it to implement 
several different cellular behaviors, including cell motility based on actin treadmilling, chemotaxis-driven cell migration, and cell proliferation in the context of cancer.

The complete report and model specifications can be found <a href="https://drive.google.com/file/d/13li3Fph4upv7blA_9rAFgvEld-MjQF-R/view" target="_blank">here</a>.

![Screenshot 1](https://cvenegasj.github.io/02712_project/cell_motility_1.png)

## Methods
The basic Cellular Potts Model (CPM) has three main characteristics:
1. It models the energy due to cell membrane contacts of neighboring cells using the adhesion matrix J, where there is a specific energy value for each pair of cell types.
2. All cells are constrained to keep an almost constant area (or volume for 3D versions)throughout state transitions.
3. Stochastic cell membrane fluctuations are simulated so that shape of cells changes while keeping their symmetry and maintaining their areas. 

Thus, in order to model the proposed cellular behaviors, we need to extend the basic version of CPM. To reproduce cell 
motility, we model the actin-elongation process close to cell membranes, and further extend this by simulating the behavior 
of cells in response to a chemoattractant. Our third extension will try to replicate the cell growth and proliferation 
phenomena observed in tumours.

## Technologies used
* The app has been written entirely in the JavaScript language (ES6 specification).
* The [p5.js](http://p5js.org/) library was used for visualization and rendering purposes.

## Demo videos
1. Actin-Driven Cell Motility CPM Simulation: <a href="https://youtu.be/GcAfP1gaiTk" target="_blank">https://youtu.be/GcAfP1gaiTk</a>
2. Chemotaxis CPM Simulation: <a href="https://youtu.be/XPplzS0niqs" target="_blank">https://youtu.be/XPplzS0niqs</a>
3. Tumor Growth CPM Simulation: <a href="https://youtu.be/PZ6Z2IzU2Js" target="_blank">https://youtu.be/PZ6Z2IzU2Js</a>
