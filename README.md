# 3GeoCombine2
![11_09_24_four_dips](https://github.com/user-attachments/assets/87e31f1f-a98d-4861-83af-62805aa22317)

The objective of this software is to analyze directional data resulting from applying a combinatorial algorithm to a displaced geological horizon. The above picture presents data stored in a regular grid. The left panels present data (all possible triangles) without elevation errors. Because there are only two elevations, our formal analysis shows that there will be identical dip directions for different triangles sharing the same edge on the flat surface. The right panels present data with added uncertainty. There is a cloud at the centre of the plot indicating the orientation of the flat surfaces. Area with lower density of observations opposite to the true dip direction can also be observed.

## Step 1 - preparing data

Prepare a text file which contains XYZ coordinates of your displaced horizon. Note that columns should be separated by a space:

1 2 4 
2.1 4 0.1
5 1 9.7

Alternatively, you can:
-download our input data from Zenodo (see reference in the paper)
-prepare synthetic data with two elevations and then add uncertainty to the elevation data (adding_error.R)

## Step 2 - applying a combinatorial algorithm

If you prepared your text file, you can now use the CPP program (main.cpp) to create all triangles based on your point data set.

## Step 3 - statistical analysis

If you have the outout from the combinatorial algorithm, you can move to statistical analysis (statistical_analysis.R) to calculate relevant directional statistics such as mean and median direction and circular dispersion.

## Appendix

Do not forget about the Zenodo repository where you can find all data used in our article!
