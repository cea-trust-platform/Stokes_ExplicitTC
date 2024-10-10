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
% RHSStokesVortexNC.m:
%
% Calcul du second membre pour resoudre le probleme de Stokes "Vortex" 
%   avec les elements finis non conformes CR ou FS 
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= F(\nu,\alpha,\beta) 
%        $F(\nu,\alpha,\beta)=-\nu(\alpha^2-1)\rho^{\alpha-2}\evec_\theta
%                               +\beta\rho^{\beta-1}\evec_\rho
%        $\rho$ est la distance a (1/2,1/2).
%    div U = 0;
%
% U(x,y)=\rho^\alpha\evec_\theta; 
% p(x,y)=\rho^\beta;
%
%
% SYNOPSIS [RHS_Ux,RHS_Uy,RHS_p,UxexNC,UyexNC,PexNC]=RHSStokesVortexNC(Ku,Mp,Bx,By,Vol)
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - NumTri(Nbtri,3) : liste de triangles (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%        - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(e,:) donne 
%                   les numeros des triangles contenant l'arete e
%		     - RefEdg(Nbedg,1) : Reference de chaque arete
%        - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2
%
% OUTPUT - RHS_Ux : second membre composante x
%        - RHS_Uy : second membre composante y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHS_Ux,RHS_Uy]=RHSStokesVortexNC()
global Nbpt CoorNeu
global Nbtri Aires NumTri RefTri
global Nbedg TriEdg
global nu alpha beta ordre
%
if (ordre==1)
  basis='CR';
end
%
if (ordre==2)
  basis='FS';
end
%
NdofU=Nbedg;
if (ordre==2)
  NP2=Nbedg+Nbpt;
  NdofU=NP2+Nbtri;
end
RHS_Ux=zeros(NdofU,1); RHS_Uy=zeros(NdofU,1);
%
ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
cnu=-nu*(alpha^2-1);
%
% Triangles qui touchent le centre
[RHS_Ux,RHS_Uy]=RHSStokesVortexPolTri(nu,alpha,beta,basis,NdofU);
%
% Triangles qui ne touchent pas le centre
TriE=find(RefTri==0); NbtriE=size(TriE,1);
for tri=1:NbtriE
  t=TriE(tri,1);
  IGLO=NumTri(t,:);
  AGLO=TriEdg(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  % volume of the element
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  X0=xyp(:,1)-0.5; Y0=xyp(:,2)-0.5;
  theta=atan2(Y0,X0);
  rho=sqrt(X0.^2+Y0.^2);
  cth=cos(theta); sth=sin(theta);
  rhoA=rho.^(alpha-2); rhoB=rho.^(beta-1);
  FxLOC = awp.*(-cnu*rhoA.*sth+beta*rhoB.*cth);
  FyLOC = awp.*( cnu*rhoA.*cth+beta*rhoB.*sth);
  % EF de CR P1NC-P0 (ordre=1)
  if (ordre==1)
    phiNC=1-2*lambda;
    UGLO=AGLO;
  end
  % EF de FS P2NC-P1 (ordre=2)
  if (ordre==2)
    lambda2=lambda.*lambda;
    phiS=2*lambda2-lambda;
    phiM(:,ILOC)=4*lambda(:,JLOC).*lambda(:,KLOC);
    phiT=2-3*sum(lambda2,2);
    phiNC=[phiS,phiM,phiT];
    UGLO=[IGLO,AGLO+Nbpt,NP2+t];
  end
  RHS_Ux(UGLO)+=sum(FxLOC.*phiNC,1)';
  RHS_Uy(UGLO)+=sum(FyLOC.*phiNC,1)';
  %
end
%