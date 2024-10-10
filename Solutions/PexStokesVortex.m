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
% PexStokesVortex.m:
%
% Resolution du probleme de Stokes "Vortex" dans un carre [0,1]*[0,1]
%   avec les elements finis P1-P0 ou P2-P1dg
%   Pb de Stokes dans un carre [0,1]*[0,1] 
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= f(Uex,pex)
%    div U = 0;
%
% U(x,y)=\rho^\alpha\,\evec_\theta; 
% p(x,y)=(\rho^\beta-\ul{\rho^\beta});
%
%
% SYNOPSIS [PexDG,PexLG,PexDG2h]=PexStokesVortex(Mu,Mp)
%
% GLOBAL - Nbtri            : nb de triangles
%        - CoorNeu(Nbpt,2)  : coordonnees (x, y) des sommets
%        - NumTri(Nbtri,3)  : liste de triangles (3 numeros de sommets)
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
%
% INPUT  - Mu,Mp : matrice de masse vitesse, matrice de masse pression
%
% OUTPUT - PexDG : la solution exacte aux points de discretisation LG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PexDG,PexLG,Pex2h]=PexStokesVortex(Mu,Mp)
%
global Nbpt CoorNeu
global Nbtri NumTri Aires TriEdg
global alpha beta ordre Vol solP
%
NdofU=size(Mu,1); NdofP=size(Mp,1);
%
RHS_P=zeros(NdofU,1);
%
%%%%%%%%%%%%%%
% Pression
%%%%%%%%%%%%%%
PexDG=zeros(NdofP,1); PexLG=zeros(NdofU,1); 
moyPex=0; Pex2h=0;
%
tdeb=1;
moy=0;
%
if (beta==1)
  moy=0.3826;
end
%
if (beta==2)
  moy=1/6;
end
%
if (beta==0.55)
  moy=0.577;
end
%
if (ordre==2)
  ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
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
  xp=xyp(:,1)-0.5; yp=xyp(:,2)-0.5;
  rho2=yp.^2+xp.^2;
  rhoB=rho2.^(beta/2); 
  PexLOC=rhoB-moy;
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
  RHS_P(UGLO)+=sum(awp.*PexLOC.*psiLG,1)';
end
%
PexLG=Mu\RHS_P;
PexLG=PexLG-RHS_P'*ones(NdofU,1)/Vol;
%%%%%%%%%%%%%%
% Pression
%%%%%%%%%%%%%%
PexDG=PexDG-moyPex*ones(NdofP,1)/Vol;
Pex2h=PexDG'*Mp*PexDG;
