I collect sample code for a project I have written in python for my PhD. 
The main script is located in Analysis/HNL.py, which calls associated libraries in Analysis/lib/   There are also bash scripts in bash/ which automate everything.
My intention is to show the type of code that I frequently develop and maintain for particle physics analyses.

The goal of this project was to look for events (or collisions) containing only a single photon (light particle) and electrons in the final state at the Large Hadron Collider. We impose cuts on the collisions that we've collected, namely we require a certain minimal energy for the photon and electron, angular direction, and so on.

Disclaimer: The goal is not for the code to run right out of the box (although with some work, one could).
The code is meant to run on large input data files which get deleted at the end of the analysis.


