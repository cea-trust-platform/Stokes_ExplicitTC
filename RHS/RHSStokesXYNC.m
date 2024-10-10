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
% RHSStokesXYNC.m:
%
% Calcul du second membre pour resoudre le probleme de Stokes "XY"
%   avec les elements finis non conformes CR ou FS
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p=2*[(X-X_0,Y-Y_0]
%    div U = 0;
%
% U(x,y)=(-Y,X);
% p(x,y)=X^3+Y^3-1/2;
%
%
% SYNOPSIS [RHS_Ux,RHS_Uy,RHS_P,Uxex_NC,Uyex_NC,Pex_NC,moyPex]=RHSStokesXYNC(Ku,Mp,Bx,By,Vol)
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : liste de triangles (3 numeros de sommets)
%        - CoorBary(Nbtri,2) : coordonnees (x, y) des barycentres des triangles
%        - Nbedg : nombre d'aretes
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%        - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(e,:) donne 
%                   les numeros des triangles contenant l'arete e
%		     - RefEdg(Nbedg,1) : Reference de chaque arete
%		     - NumEdg(Nbedg,2) : Numeros des sommets de chaque arete
%        - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2
%
% INPUT  - Ku,Bx,By : matrice de raideur et de couplage (relevement CLDNH)
%
% OUTPUT - RHS_Ux : second membre composante x
%        - RHS_Uy : second membre composante y
%        - RHS_P  : second membre pression (relevement CLDNH)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHS_Ux,RHS_Uy]=RHSStokesXYNC()
global Nbpt CoorNeu RefNeu
global Nbtri Aires NumTri RefTri
global Nbedg TriEdg RefEdg EdgTri NumEdg
global ordre
%
%
Ndof=Nbedg;
if (ordre==2)
  NP2=Nbedg+Nbpt;
  Ndof=NP2+Nbtri;
end
RHS_Ux=zeros(Ndof,1); RHS_Uy=zeros(Ndof,1);
%
ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
%
for t=1:Nbtri
  IGLO=NumTri(t,:);
  AGLO=TriEdg(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  % volume of the element
  aire=Aires(t);
  % Integration points
  [xyp,wp,lmbd,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  FxLOC = 3*awp.*(xyp(:,1).^2); FyLOC = 3*awp.*(xyp(:,2).^2);
  % EF de CR P1NC-P0 (ordre=1)
  if (ordre==1)
    phiNC=1-2*lmbd;
    UGLO=AGLO;
  end
  % EF de FS P2NC-P1 (ordre=2)
  if (ordre==2)
    lmbd2=lmbd.*lmbd;
    phiS=2*lmbd2-lmbd;
    phiM(:,ILOC)=4*lmbd(:,JLOC).*lmbd(:,KLOC);
    phiT=2-3*sum(lmbd2,2);
    phiNC=[phiS,phiM,phiT];
    UGLO=[IGLO,AGLO+Nbpt,NP2+t];
  end
  RHS_Ux(UGLO)+=sum(FxLOC.*phiNC,1)';
  RHS_Uy(UGLO)+=sum(FyLOC.*phiNC,1)';
end
%