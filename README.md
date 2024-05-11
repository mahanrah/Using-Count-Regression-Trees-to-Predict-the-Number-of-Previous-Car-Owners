# Using Count Regression Trees to Predict the Number of Previous Car Owners

ABSTRACT

The purpose of this project was to research count regression trees, and to apply them to a dataset of Irish used cars. The aim of this application was to predict the number of previous owners based on variables such as vehicle registration year and price. An in depth understanding of the CORE method and its applications was developed through a review of existing literature. The chosen dataset was then analysed, imputed and prepared for inclusion in the model. The  Analysis found that the model was biased toward lower numbers of previous owners, and did not perform well for higher values of this variable.

The R file in this repository contains functions to create a count regression tree for the Irish_UsedCars dataset, which can be used to predict the number of previous owners of a car.

The CORE method is taken from Liu, Lin, and Shih (2020), which is accessible from this link: https://link.springer.com/article/10.1007/s11634-019-00358-7

INSTRUCTIONS

Download the Irish_UsedCars dataset from Kaggle (https://www.kaggle.com/datasets/ravi23231/irish-usedcars) or any other count dataset, download the required packages, and run the functions

REQUIRED PACKAGES

MASS, pscl, glmmTMB

LICENCE

This project is licensed under the MIT License - see the  file for details.

AUTHORS

Áine Sinnott, Buse Suer, Cera Kenny, Dylan Mangan, Laura Kealy, Matthew Hanrahan, Phoebe Morgan, Will Robbins

CONTACT

mahanrah@tcd.ie
