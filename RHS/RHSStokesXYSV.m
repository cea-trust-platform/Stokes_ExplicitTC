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
% RHSStokesXYSV.m:
%
% Calcul du second membre pour resoudre le probleme de Stokes "U=(-Y,X)" dans un carre [0,1]*[0,1]
%   avec les elements finis P1-P0 ou P2-P1dg
%   -nu*Delta U + grad p=2*[(X-X_0,Y-Y_0]
%    div U = 0;
%
% U(x,y)=(-Y,X);
% p(x,y)=X^3+Y^3-1/2;
%
% SYNOPSIS [RHS_Ux,RHS_Uy,RHS_p]=RHSStokesXYSV(Ku,Mp)
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : liste de triangles (3 numeros de sommets)
%        - CoorBary(Nbtri,2) : coordonnees (x, y) des barycentres des triangles
%        - Nbedg : nombre d'aretes
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - RefEdg(Nbedg,1) : Reference de chaque arete
%        - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2
%
% OUTPUT - RHS_Ux : second membre vitesse composante x
%        - RHS_Uy : second membre vitesse composante y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHS_Ux,RHS_Uy]=RHSStokesXYSV()
global Nbpt CoorNeu
global Nbtri Aires NumTri
global Nbedg TriEdg
global ordre
%
% Second membre 
% conservation de la quantite de mouvement
%
Ndof=Nbpt;
if (ordre==2)
  Ndof+=Nbedg;
end
%
RHS_Ux=zeros(Ndof,1); RHS_Uy=zeros(Ndof,1);
%
ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
for t=1:Nbtri
  IGLO=NumTri(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  AGLO=TriEdg(t,:);
  % volume of the element
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  FxLOC = 3*awp.*(xyp(:,1).^2); FyLOC = 3*awp.*(xyp(:,2).^2);
  %
  % EF P1-P0 (ordre=1)
  if (ordre==1)
    phiLG=lambda;
    UGLO=IGLO;
  end
  % EF P2-P1dg (ordre=2)
  if (ordre==2)
    lambda2=lambda.*lambda;
    phiS=2*lambda2-lambda;
    phiM(:,ILOC)=4*lambda(:,JLOC).*lambda(:,KLOC);
    phiLG=[phiS,phiM];
    UGLO=[IGLO,AGLO+Nbpt];
  end
  RHS_Ux(UGLO)+=sum(FxLOC.*phiLG,1)';
  RHS_Uy(UGLO)+=sum(FyLOC.*phiLG,1)';
end
%