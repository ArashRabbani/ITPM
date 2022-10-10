% In this script you are comparing the absolute permeabilities of a
% tube/throat with arbitrary cross-section which is loaded from an image in
% the current directory

clear; clc; close all;
A=imread('B.png'); % Input image name, yuo can change this 
if ndims(A)==3; A=rgb2gray(A); end; A=im2bw(double(A),graythresh(A)); % converting the input image into binary 
A=imresize(A,.07,'nearest'); % Make the input image smaller with desired ratio

Plot=1;  % if set to 1 it will show the LBM convergence plots and else it will be hidden
Res=1;   % This is the spatial correlation of the images 

Method='EMP';   % Use Empirical Correlation
AbsolutePerm1=ITPM(A,Res,Method,Plot);

Method='LBM';   % Use Lattice Boltzmann simulation
AbsolutePerm2=ITPM(A,Res,Method,Plot);

Method='ANN1P'; % Use Neural Network model with 1 input paramters
AbsolutePerm3=ITPM(A,Res,Method,Plot);

Method='ANN7P'; % Use Neural Network model with 7 input paramters
AbsolutePerm4=ITPM(A,Res,Method,Plot);

