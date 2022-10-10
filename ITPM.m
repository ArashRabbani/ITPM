function [AbsolutePerm]=ITPM(A,Res,Method,Plot)
% Image-based tube/throat permeability model is a mean to find the absolute
% permeability of tube with arbitrary cross-section this function can use 4
% mthods for estimating the absolute permeability: 1) Latice Boltzmann
% simulation, 2) An artificial neural network with 1 input paramter , 3)
% Another artificial neural network with 7 input paramter and , 4) an
% empirical correlation which uses the average distance values of the
% transformed input images

% Inputs:
% A: is a binary image in which void space is 0 and solid space is 1, this
% image shows the cross-section of the throat/tube
% Res: is the spatial resolution and it is expressed as micron/pixel
% Method: asks that what method you wanted to use for permeability
% calculation the values could be : LBM, EMP, ANN1P, and ANN7P.
% Plot: when put as 1 it will shows the LBM convergence charts and if set
% to zero it wont

% Output: Absolute Permeability of throat/tube in Darcy
% The LBM section is adopted from this source: 
% Haslam, I. W., Crouch, R. S., & Seaïd, M. (2008). Coupled finite element–lattice
% Boltzmann analysis. Computer Methods in Applied Mechanics and Engineering, 197(51-52), 4505-4511.

% Note: You will need neural fitting and image processing toolboxes to use
% this code

if strcmp(Method,'LBM')
for j=1:6; GEO(j,:,:)=A; end 
nx=size(GEO,1); ny=size(GEO,2); nz=size(GEO,3);
omega=1; density=1 ;t1=1/3; t2=1/18; t3=1/36;
F=repmat(density/19,[nx ny nz 19]); FEQ=F; matsize=nx*ny*nz;
CI=[0:matsize:matsize*19];
ON=find(GEO); %matrix offset of each Occupied Node
TO_REFLECT=[ON+CI(2) ON+CI(3) ON+CI(4) ON+CI(5)	ON+CI(6) ON+CI(7) ON+CI(8) ...
 ON+CI(9) ON+CI(10) ON+CI(11) ON+CI(12) ON+CI(13) ON+CI(14) ON+CI(15) ...
 ON+CI(16) ON+CI(17) ON+CI(18) ON+CI(19)];
REFLECTED=[ON+CI(3) ON+CI(2) ON+CI(5) ON+CI(4) ON+CI(7) ON+CI(6) ON+CI(11) ...
 ON+CI(10) ON+CI(9) ON+CI(8) ON+CI(15) ON+CI(14) ON+CI(13) ON+CI(12) ...
 ON+CI(19) ON+CI(18) ON+CI(17) ON+CI(16)];
avu=1; prevavu=1; ts=0; deltaU=1e-7; numactivenodes=sum(sum(sum(1-GEO)));
if Plot==1
figure('units','normalized','outerposition',[0 0 1 1]);  subplot(2,2,1); hold on; xlabel('Iterations'); ylabel('Logarithm of Permeability Relative Error');
subplot(2,2,3); hold on; xlabel('Iterations'); ylabel('Aboslute Permeability (Darcy)'); 
end
Error=1; i=1; 
% while (ts<4000 & 1e-6<abs((prevavu-avu)/avu)) | ts<100
% while (ts<4000 && 1e-6<PermError(i) || ts<100)
while (ts<4000 && 1e-6<Error) || ts<100
 % Propagate
 %nearest-neighbours
 F(:,:,:,2)=F(:,:,[nz 1:nz-1],2);
 F(:,:,:,3)=F(:,:,[2:nz 1],3);
 F(:,:,:,4)=F(:,[ny 1:ny-1],:,4);
 F(:,:,:,5)=F(:,[2:ny 1],:,5);
 F(:,:,:,6)=F([nx 1:nx-1],:,:,6);
 F(:,:,:,7)=F([2:nx 1],:,:,7);
 %next-nearest neighbours
 F(:,:,:,8)= F([nx 1:nx-1],[ny 1:ny-1],:,8);
 F(:,:,:,9)= F([nx 1:nx-1],[2:ny 1],:,9);
 F(:,:,:,10)=F([2:nx 1],[ny 1:ny-1],:,10);
 F(:,:,:,11)=F([2:nx 1],[2:ny 1],:,11);
 F(:,:,:,12)=F([nx 1:nx-1],:,[nz 1:nz-1],12);
 F(:,:,:,13)=F([nx 1:nx-1],:,[2:nz 1],13);
 F(:,:,:,14)=F([2:nx 1],:,[nz 1:nz-1],14);
 F(:,:,:,15)=F([2:nx 1],:,[2:nz 1],15);
 F(:,:,:,16)=F(:,[ny 1:ny-1],[nz 1:nz-1],16);
 F(:,:,:,17)=F(:,[ny 1:ny-1],[2:nz 1],17);
 F(:,:,:,18)=F(:,[2:ny 1],[nz 1:nz-1],18);
 F(:,:,:,19)=F(:,[2:ny 1],[2:nz 1],19);
 BOUNCEDBACK=F(TO_REFLECT); %Densities bouncing back at next timestep
 % Relax; calculate equilibrium state (FEQ) with equivalent speed and density to F
 DEN = sum(F,4);
 UX=(sum(F(:,:,:,[6 8 9 12 13]),4)-sum(F(:,:,:,[7 10 11 14 15]),4))./DEN;
 UY=(sum(F(:,:,:,[4 8 10 16 17]),4)-sum(F(:,:,:,[5 9 11 18 19]),4))./DEN;
 UZ=(sum(F(:,:,:,[2 12 14 16 18]),4)-sum(F(:,:,:,[3 13 15 17 19]),4))./DEN;
 UX(1,:,:)=UX(1,:,:)+deltaU; %Increase inlet pressure
 UX(ON)=0; UY(ON)=0; UZ(ON)=0; DEN(ON)=0; U_SQU=UX.^2+UY.^2+UZ.^2;
 U8=UX+UY;U9=UX-UY;U10=-UX+UY;U11=-U8;U12=UX+UZ;U13=UX-UZ;
 U14=-U13;U15=-U12;U16=UY+UZ;U17=UY-UZ;U18=-U17;U19=-U16;
 % Calculate equilibrium distribution: stationary
 FEQ(:,:,:,1)=t1*DEN.*(1-3*U_SQU/2);
 % nearest-neighbours
 FEQ(:,:,:,2)=t2*DEN.*(1 + 3*UZ + 9/2*UZ.^2 - 3/2*U_SQU);
 FEQ(:,:,:,3)=t2*DEN.*(1 - 3*UZ + 9/2*UZ.^2 - 3/2*U_SQU);
 FEQ(:,:,:,4)=t2*DEN.*(1 + 3*UY + 9/2*UY.^2 - 3/2*U_SQU);
 FEQ(:,:,:,5)=t2*DEN.*(1 - 3*UY + 9/2*UY.^2 - 3/2*U_SQU);
 FEQ(:,:,:,6)=t2*DEN.*(1 + 3*UX + 9/2*UX.^2 - 3/2*U_SQU);
 FEQ(:,:,:,7)=t2*DEN.*(1 - 3*UX + 9/2*UX.^2 - 3/2*U_SQU);
 % next-nearest neighbours
 FEQ(:,:,:,8) =t3*DEN.*(1 + 3*U8 + 9/2*(U8).^2 - 3*U_SQU/2);
 FEQ(:,:,:,9) =t3*DEN.*(1 + 3*U9 + 9/2*(U9).^2 - 3*U_SQU/2);
 FEQ(:,:,:,10)=t3*DEN.*(1 + 3*U10 + 9/2*(U10).^2 - 3*U_SQU/2);
 FEQ(:,:,:,11)=t3*DEN.*(1 + 3*U11 + 9/2*(U11).^2 - 3*U_SQU/2);
 FEQ(:,:,:,12)=t3*DEN.*(1 + 3*U12 + 9/2*(U12).^2 - 3*U_SQU/2);
 FEQ(:,:,:,13)=t3*DEN.*(1 + 3*U13 + 9/2*(U13).^2 - 3*U_SQU/2);
 FEQ(:,:,:,14)=t3*DEN.*(1 + 3*U14 + 9/2*(U14).^2 - 3*U_SQU/2);
 FEQ(:,:,:,15)=t3*DEN.*(1 + 3*U15 + 9/2*(U15).^2 - 3*U_SQU/2);
 FEQ(:,:,:,16)=t3*DEN.*(1 + 3*U16 + 9/2*(U16).^2 - 3*U_SQU/2);
 FEQ(:,:,:,17)=t3*DEN.*(1 + 3*U17 + 9/2*(U17).^2 - 3*U_SQU/2);
 FEQ(:,:,:,18)=t3*DEN.*(1 + 3*U18 + 9/2*(U18).^2 - 3*U_SQU/2);
 FEQ(:,:,:,19)=t3*DEN.*(1 + 3*U19 + 9/2*(U19).^2 - 3*U_SQU/2);
 F=omega*FEQ+(1-omega)*F;
 F(REFLECTED)=BOUNCEDBACK;
 prevavu=avu;avu=sum(sum(sum(UX)))/numactivenodes; ts=ts+1;

  
 % Permeability Calculation
 rho = sum(F,4);
 temp=UX(end,:,:); temp=temp(:); temp(temp==0)=[]; Ubar=mean(temp);
 Ubar = mean(UX(:));
 P = rho/3;
 Pin = P(end-1,:,:);
 Pout = P(end,:,:);
 dP = mean(Pout(:)) - mean(Pin(:));
 viscosity = (1/omega - 1/2)/3;
 K = - Ubar * viscosity / dP;
 K=K*Res^2/.9869; %Darcy
 Perm(i)=K; 
 if i>1; Error= abs(Perm(i)- Perm(i-1))/Perm(i);  PermError(i)=Error;  end
 
 if Plot==1 && i>10 && mod(i,10)==1;
 subplot(2,2,1); plot([i-10,i],log10([PermError(i-10),PermError(i)]),'k'); title(['Permeability Relative Error =' num2str(Error)]);
 subplot(2,2,3); plot([i-10,i],([Perm(i-10),Perm(i)]),'k'); title(['Absolute Permeability =' num2str(round(K,3)) ' (Darcy)']);
 
subplot(2,2,[2,4]); 
Uout=squeeze(UX(end,:,:)); 
surf((Uout).*1e6);
c = colorbar; colormap jet;
c.Label.String = 'Velocity (\mum/s)'; axis equal square tight off;
h.Children(1).FontSize=14; h.Children(1).FontWeight='bold'; title('Velocity Distribution');
  drawnow; end
   i=i+1;
end
AbsolutePerm=Perm(end);
end
if strcmp(Method,'EMP')
  B=bwdist(A); B=B(:); B(B==0)=[];   
  Di=mean(B); AbsolutePerm=Di^2*1.3042-1.0851*Di-0.1;
  AbsolutePerm=AbsolutePerm*Res^2/.9869; %Darcy
end
if strcmp(Method,'ANN1P')
    BW=1-A;
    SE2=zeros(3,3); SE2(2,2)=1; SE2(2,3)=1; SE2(2,1)=1; SE2(3,2)=1; SE2(1,2)=1;
    reg=regionprops(BW,'ConvexArea','Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength','Solidity');
    Area=sum(sum(BW));
    Perimeter=sum(sum(bwperim(imdilate(BW,SE2),4)));
    Axisratio=reg.MajorAxisLength/reg.MinorAxisLength;
    EqDiam=reg.EquivDiameter;
    Solidity=reg.Solidity;
    Ecc=reg.Eccentricity;
    B=bwdist(A); B=B(:); B(B==0)=[];   
    MeanDist=mean(B);
    load('ANN1P.mat')
    AbsolutePerm=net(MeanDist)*Res^2/.9869; %Darcy
end
if strcmp(Method,'ANN7P')
    BW=1-A;
    SE2=zeros(3,3); SE2(2,2)=1; SE2(2,3)=1; SE2(2,1)=1; SE2(3,2)=1; SE2(1,2)=1;
    reg=regionprops(BW,'ConvexArea','Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength','Solidity');
    Area=sum(sum(BW));
    Perimeter=sum(sum(bwperim(imdilate(BW,SE2),4)));
    Axisratio=reg.MajorAxisLength/reg.MinorAxisLength;
    EqDiam=reg.EquivDiameter;
    Solidity=reg.Solidity;
    Hydrad=2*Area/Perimeter;
    B=bwdist(A); B=B(:); B(B==0)=[];   
    MeanDist=mean(B);
    load('ANN7P.mat')
    AbsolutePerm=net([Area,Perimeter,Axisratio,EqDiam,Solidity,Hydrad,MeanDist]')*Res^2/.9869; %Darcy
end
end
