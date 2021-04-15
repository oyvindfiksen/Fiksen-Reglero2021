# Fiksen-Reglero2021
Atlantic bluefin tuna larvae growth and fitness

This repository contains all files needed to run the model and generate the figures presented in the paper (Fiksen & Reglero: "Atlantic bluefin tuna spawn early to avoid metabolic meltdown in larvae") and the Appendix. 
The main source code 'FiguresforPaper.py' is written in Python 7.3 - by Øyvind Fiksen (oyvind.fiksen@uib.no).

The code calculates the fitness (survival chance) of egg of Atlantic Bluefin Tuna (ABFT) spawned at each day of the year, from day of birth to the flexion stage. 

The code also generates a number of plots and .svg - files that together form the basis of figures in the manuscript. The model is explained in detail in the Supplementary. 

Two text-files provide input data on daylength and daily average temperture to the model (HoursofLight.txt and AverageTempData_NOAA.txt), and these are included here. Other published data are embedded into the code - the data on larval densities from field surveys (from Reglero & al 2018)

A number of .npy files provide model results from a range of model runs with different parameter or settings of food and temperature. These are needed to generate the figures in the paper and the appendix.

Øyvind Fiksen, April 2021
Department of Biological Sciences, University of Bergen, Norway
http://bio.uib.no/te/of/
oyvind.fiksen@uib.no
