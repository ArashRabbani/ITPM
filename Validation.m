% This script compares the claculated absolute permeabilities using 4
% methods of ITPM with the analytical answer for a cylindrical tube 

clear; clc; close all;
CircleRad=7; % pixels % you can change this radius
A=1-sph(CircleRad); 

Plot=0;  % if set to 1 it will show the LBM convergence plots and else it will be hidden
Res=1;   % This is the spatial correlation of the images

Method='EMP';   % Use Empirical Correlation
AbsolutePerm1=ITPM(A,Res,Method,Plot); %Darcy

Method='LBM';   % Use Lattice Boltzmann simulation
AbsolutePerm2=ITPM(A,Res,Method,Plot); %Darcy

Method='ANN1P'; % Use Neural Network model with 1 input paramter
AbsolutePerm3=ITPM(A,Res,Method,Plot); %Darcy

Method='ANN7P'; % Use Neural Network model with 7 input paramters
AbsolutePerm4=ITPM(A,Res,Method,Plot); %Darcy

% Calculate the absolute permeability of a tube based on the Hagen Poiseuille equation K=r^2/8   
AbsolutePerm5=(CircleRad+.25)^2/8*Res^2/.9869; %Darcy

figure; 
bar([AbsolutePerm1,AbsolutePerm2,AbsolutePerm3,AbsolutePerm4,AbsolutePerm5]);
ylabel('Absolute Permeability (Darcy)'); set(gca,'xticklabel',{'Empirical','LBM','ANN1P','ANN7P','Analytical'})