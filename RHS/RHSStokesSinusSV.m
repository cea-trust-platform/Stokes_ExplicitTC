%{
/****************************************************************************
* Copyright (c) 2022, CEA
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
% RHSStokesSinusSV.m:
%
% Resolution du probleme de Stokes "Sinus" dans un carre [0,1]*[0,1]
%   avec les elements finis P1-P0 ou P2-P1dg
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= npi*[sin(npi*y)*{(bp-2*nu*npi)cos(npi*x)+nu*npi} ; 
%                              sin(npi*x)*{(bp+2*nu*npi)cos(npi*y)-nu*npi)}]
%    div U = 0;
%
% U(x,y)=[(1-cos(npi*x))*sin(npi*y) ; 
%         (cos(npi*y)-1)*sin(npi*x)]; 
% p(x,y)=bp*sin(npi*x)*sin(npi*y);
%
% Calcul du second membre
%
% SYNOPSIS [RHSx,RHSy]=RHSStokesSinusSV()
%
% GLOBAL - Nbpt : nombre de sommets
%        - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - Nbtri : nombre de triangles
%        - Aires : aire des triangles
%        - ordre : 1 ou 2
%
% OUTPUT - RHS_Ux : second membre vitesse composante x
%        - RHS_Uy : second membre vitesse composante y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHS_Ux,RHS_Uy]=RHSStokesSinusSV( )
global Nbpt CoorNeu
global Nbtri Aires NumTri
global Nbedg TriEdg
global ordre
global nu npi npi2 bp
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
npi2nu=npi2*nu;
cx=npi*bp-2*npi2nu;
cy=npi*bp+2*npi2nu;
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
  sinxy=sin(npi*xyp); cosxy=cos(npi*xyp);
  FxLOC = awp.*sinxy(:,2).*(cx*cosxy(:,1)+npi2nu);
  FyLOC = awp.*sinxy(:,1).*(cy*cosxy(:,2)-npi2nu);
  %
  % EF P1-P0 (ordre=1)
  if (ordre==1)
    phiLG=lambda; UGLO=IGLO;

  end
  % EF P2-P1dg (ordre=2)
  if (ordre==2)
    lambda2=lambda.*lambda;
    phiS=2*lambda2-lambda;
    phiM(:,ILOC)=4*lambda(:,JLOC).*lambda(:,KLOC);
    phiLG=[phiS,phiM]; UGLO=[IGLO,AGLO+Nbpt];
  end
  RHS_Ux(UGLO)+=sum(FxLOC.*phiLG,1)';
  RHS_Uy(UGLO)+=sum(FyLOC.*phiLG,1)';
end
%