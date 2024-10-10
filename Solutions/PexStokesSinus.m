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
% PexStokesSinus.m:
%
% Resolution du probleme de Stokes "Sinus" dans un carre [0,1]*[0,1] 
%   avec les elements finis P1-P0 ou P2-P1dg
%   -nu*Delta U + grad p= npi*[sin(npi*y)*{(bp-2*nu*npi)cos(npi*x)+nu*npi} ; 
%                              sin(npi*x)*{(bp+2*nu*npi)cos(npi*y)-nu*npi)}]
%    div U = 0;
%
% U(x,y)=[(1-cos(npi*x))*sin(npi*y) ; 
%         (cos(npi*y)-1)*sin(npi*x)]; 
% p(x,y)=bp*sin(npi*x)*sin(npi*y);
%
% Calcul de la solution
%
% SYNOPSIS [PexDG,moyPexDG,PexLG,Pex2h]=PexStokesSinus(Mu,Mp)
%
% GLOBAL - Nbtri            : nb de triangles
%        - CoorNeu(Nbpt,2)  : coordonnees (x, y) des sommets
%        - NumTri(Nbtri,3)  : liste de triangles (3 numeros de sommets)
%        - Nbedg            : nb d'aretes
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
%
% INPUT  - Mu,Mp : matrice de masse vitesse, matrice de masse presison
%
% OUTPUT - PexDG : la solution exacte aux points de discretisation NC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PexDG,moyPexDG,PexLG]=PexStokesSinus(Mu,Mp)
global Nbpt CoorNeu
global Nbtri NumTri Aires TriEdg
global ordre npi bp
global Vol
%
NdofU=size(Mu,1); NdofP=size(Mp,1);
%
[Mass1,invMass1,Mass01,Mass2,invMass2,Mass02,Mass12]=MassLGTriangle();
%
RHS_Uxex=zeros(NdofU,1); RHS_Uyex=RHS_Uxex; RHS_PexLG=RHS_Uxex;
RHS_Ux=RHS_Uxex; RHS_Uy=RHS_Uyex;
%
ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
%
%%%%%%%%%%%%%%
% Pression
%%%%%%%%%%%%%%
PexDG=zeros(NdofP,1); PexLG=zeros(NdofU,1); 
moyPexDG=0;
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
  npixyp=npi*xyp;
  cosxy=cos(npixyp); sinxy=sin(npixyp);
  PexLOC=bp*sinxy(:,1).*sinxy(:,2);
  PexT=wp*PexLOC;
  moyPexDG+=aire*PexT;
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
  RHS_PexLG(UGLO)+=sum(awp.*PexLOC.*psiLG,1)';
end
%
%%%%%%%%%%%%%%
% Pression
%%%%%%%%%%%%%%
PexDG=PexDG-moyPexDG*ones(NdofP,1)/Vol;
%
PexLG=Mu\RHS_PexLG;
moyPexLG=ones(1,NdofU)*RHS_PexLG;
PexLG=PexLG-moyPexLG*ones(NdofU,1)/Vol;
