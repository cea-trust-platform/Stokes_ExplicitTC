%{
/****************************************************************************
* Copyright (c) 2023, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot CEA
%
% PiRT1_abdXY.m:
%
% Calcul du second membre pour resoudre le probleme de Stokes
%   avec les elements finis non conformes FS 
%   et projection sur des fonctions-test du second membre sur les EF RT1.
%
% SYNOPSIS [abdX,abdY,matB]=PiRT1_abdXY(aire,CoorNeuT,EdgNormT)
%
% Pour tester aire=0.5; EdgNormT=[1,1;-1,0;0,-1];CoorNeuT=[0,0;1,0;0,1];     
%   
% INPUT - aire : aire du triangle
%       - CoorNeuT : coordonnees locales des sommets
%       - EdgNormT : coordonnees locales des vecteurs face-normale
%
% OUTPUT - [abdX,abdY,matB] : matrices necessaires pour calculer les projections des fonctions de base
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [abdX,abdY,matB]=PiRT1_abdXY(aire,CoorNeuT,EdgNormT)
% Matrice des integrales de $\Pi_{RT1}$
inv3=1/3; inv6=1/6; inv12=1/12; inv24=1/24;
M=[inv3  inv6  inv3 inv6  inv3  inv6  0.5  0.5;
   inv6  inv3  inv6 inv3  inv6  inv3  0.5  0.5;
      0 -inv3     0    0     0     0 -0.5    0;
      0 -inv6     0    0     0     0 -0.5    0;
      0     0 -inv6    0     0     0    0 -0.5;
      0     0 -inv3    0     0     0    0 -0.5;
   inv6  inv6     0    0 inv12 inv24  0.5    0;
      0     0  inv6 inv6 inv24 inv12    0  0.5];
%
invM=[ -4 -2 -2 12  6 -4  16   8 ;
        0  0 -6  6  0  0   0   0 ;
        0  0  0  0  6 -6   0   0 ;
       -2 -4 -4  6 12 -2   8  16 ;
        8  0  0 -8 -8  8 -16  -8 ;
        0  8  8 -8 -8  0  -8 -16 ;
        0  0  2 -4  0  0   0   0 ;
        0  0  0  0 -4  2   0   0 ];
% Matrice des integrales des fonctions de bases
Phi=zeros(8,7);
Phi(1,2)=inv6; Phi(1,4)=inv3; % F0, x
Phi(2,3)=inv6; Phi(2,4)=inv3; % F0, y
Phi(3,3)=inv6; Phi(3,5)=inv3; % F1, y
Phi(4,1)=inv6; Phi(4,5)=inv3; % F1, 1-y
Phi(5,1)=inv6; Phi(5,6)=inv3; % F2, 1-x
Phi(6,2)=inv6; Phi(6,6)=inv3; % F2, x
 % K, ex et ey
Phi(7:8,4:6)=inv6*ones(2,3); 
Phi(7:8,7)=0.25*ones(2,1);
%
% Matrice de la transformation affine et de son inverse
matB=[CoorNeuT(2,:)-CoorNeuT(1,:);CoorNeuT(3,:)-CoorNeuT(1,:)]';
inv2A=1/(2*aire);
matC=-inv2A*[EdgNormT(2,:);EdgNormT(3,:)];
Cxy=zeros(8,2);
for i=1:3
  ii=2*i;
  Cxy(ii,:)=inv2A*EdgNormT(i,:); Cxy(ii-1,:)=Cxy(ii,:);
end
Cxy(7:8,1:2)=matC;
% boucle sur les fonctions de base
abdX=invM*(Cxy(:,1).*Phi);
abdY=invM*(Cxy(:,2).*Phi);
