%{
/****************************************************************************
* Copyright (c) 2024, CEA
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
% PexStokesXY.m:
%
% Resolution du probleme de Stokes "U=(-Y,X)" dans un carre [0,1]*[0,1]
%   avec les elements finis P1-P0 ou P2-P1dg
%   -nu*Delta U + grad p=2*[(X-X_0,Y-Y_0]
%    div U = 0;
%
% U(x,y)=(-(Y-Y0),X-X_0);
% p(x,y)=(X-X_0)^2+(Y-Y_0)^2-1/6;
%
% Calcul de la solution
%
% SYNOPSIS [PexDG,moyPexDG,PexLG,PexDG2h]=PexStokesXY(Mu,Mp)
%
% GLOBAL - Nbtri            : nb de triangles
%        - CoorNeu(Nbpt,2)  : coordonnees (x, y) des sommets
%        - CoorMil(Nbedg,2) : coordonnees (x, y) des milieux d'aretes
%        - CoorBary(Nbtri,2): coordonnees (x, y) des barycentres des triangles
%        - NumTri(Nbtri,3)  : liste de triangles (3 numeros de sommets)
%        - Nbedg            : nb d'aretes
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
%
% INPUT  - Mu,Mp : matrice de masse et de raideur vitesse, matrice de masse presison
%
% OUTPUT - PexDG,moyPexDG,PexLG,Pex2h : la solution exacte aux points de discretisation DG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PexDG,moyPexDG,PexLG,Pex2h]=PexStokesXY(Mu,Mp)
%
global Nbpt CoorNeu
global Nbtri CoorBary NumTri Aires TriEdg
global Nbedg CoorMil
global ordre
global Vol
%
NdofU=size(Mu,1); NdofP=size(Mp,1);
PexDG=zeros(NdofP,1);
%
%%%%%%%%%%%%%%
% Pression
%%%%%%%%%%%%%%
%
moyPexDG=0; RHS_p=zeros(NdofU,1);
%
ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
%
moy=1/2; moyPex=0;
%
if (ordre==2)
  tdeb=1;
  invMass1 = [9 -3 -3;-3 9 -3;-3 -3 9];
end
%
for t=1:Nbtri
  IGLO=NumTri(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  PexLOC=xyp(:,1).^3+xyp(:,2).^3-moy;
  PexT=wp*PexLOC;
  moyPex+=aire*PexT;
  if (ordre==1)
    psiLG=lambda; UGLO=IGLO;
    PexDG(t,1)=PexT;
  end
  %
  if (ordre==2)
    AGLO=TriEdg(t,:);
    tfin=tdeb+2;
    lambda2=lambda.*lambda;
    psiS=2*lambda2-lambda;
    psiM(:,ILOC)=4*lambda(:,JLOC).*lambda(:,KLOC);
    psiLG=[psiS,psiM]; UGLO=[IGLO,AGLO+Nbpt]; 
    PexDG(tdeb:tfin,1)=invMass1*sum(wp'.*PexLOC.*lambda,1)';
    tdeb=tfin+1;
  endif
  RHS_p(UGLO)+=sum(awp.*PexLOC.*psiLG,1)';
end
%
PexLG=Mu\RHS_p;
PexLG=PexLG-RHS_p'*ones(NdofU,1)/Vol;
%%%%%%%%%%%%%%
% Pression
%%%%%%%%%%%%%%
PexDG=PexDG-moyPex*ones(NdofP,1)/Vol;
Pex2h=PexDG'*Mp*PexDG;