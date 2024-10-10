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
% RHSStokesVortexSV.m:
%
% Calcul du second membre pour resoudre le probleme de Stokes "Vortex" dans un carre [0,1]*[0,1]
%   avec les elements finis P1-P0 ou P2-P1dg
%   -nu*Delta U + grad p= F(\nu,\alpha,\beta) 
%        $F(\nu,\alpha,\beta)=-\nu(\alpha^2-1)\rho^{\alpha-2}\evec_\theta
%                               +\beta\rho^{\beta-1}\evec_\rho
%        $\rho$ est la distance a (1/2,1/2).
%    div U = 0;
%
% U(x,y)=\rho^\alpha\evec_\theta; 
% p(x,y)=\rho^\beta;
%
% SYNOPSIS [RHS_Ux,RHS_Uy,RHS_p]=RHSStokesVortexSV(Ku,Mp,Bx,By,Pex,moyPex,Uxex,Uyex)
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
% OUTPUT - RHS_Ux : second membre composante x
%        - RHS_Uy : second membre composante y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHS_Ux,RHS_Uy]=RHSStokesVortexSV()
global Nbpt CoorNeu RefNeu
global Nbtri Aires NumTri CoorBary RefTri SomOpp
global Nbedg TriEdg RefEdg EdgTri NumEdg LgEdg2 CoorMil
global nu alpha beta ordre
%
basis='LG';
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
% cnu=0 si alpha=1
cnu=nu*(alpha^2-1);
%
% Triangles qui touchent le centre
[RHS_Ux,RHS_Uy]=RHSStokesVortexPolTri(nu,alpha,beta,basis,Ndof);
%
% Triangles qui ne touchent pas le centre
TriE=find(RefTri==0); NbtriE=size(TriE,1);
for tri=1:NbtriE
  t=TriE(tri,1);
%for t=1:Nbtri
  IGLO=NumTri(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  AGLO=TriEdg(t,:);
  % volume of the element
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  X0=xyp(:,1)-0.5; Y0=xyp(:,2)-0.5;
  th0=atan2(Y0,X0);
  rho=sqrt(X0.^2+Y0.^2);
  cth=cos(th0); sth=sin(th0);
  DeltaUxLoc=0; DeltaUyLoc=0;
  if (alpha~=1)
    rhoA=rho.^(alpha-2);
    DeltaUxLoc= cnu*rhoA.*sth;
    DeltaUyLoc=-cnu*rhoA.*cth;
  end
  if (beta==2)
    GpxLoc=2*X0; GpyLoc=2*Y0;
  else
    rhoB=rho.^(beta-1);
    GpxLoc=beta*rhoB.*cth;
    GpyLoc=beta*rhoB.*sth;
  end
  FxLOC = awp.*(DeltaUxLoc+GpxLoc);
  FyLOC = awp.*(DeltaUyLoc+GpyLoc);
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