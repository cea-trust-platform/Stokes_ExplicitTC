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
% SolutionStokesSinusSV.m:
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
% SYNOPSIS [Uxex,Uyex,Pex,PexLG,Uex2h,GUex2h,Pex2h]=SolutionStokesSinusSV(Mu,Ku,DofUB)
%
% GLOBAL - Nbtri            : nb de triangles
%        - CoorNeu(Nbpt,2)  : coordonnees (x, y) des sommets
%        - NumTri(Nbtri,3)  : liste de triangles (3 numeros de sommets)
%        - Nbedg            : nb d'aretes
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
%
% INPUT  - Mu,Ku : matrice de masse et de raideur vitesse
%
% OUTPUT - Uxex,Uyex : la solution exacte aux points de discretisation LG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxex,Uyex]=SolutionStokesSinusSV(Mu,Ku)
global Nbpt CoorNeu
global Nbtri NumTri Aires TriEdg
global ordre npi bp
global Vol
%
NdofU=size(Mu,1);
%
[Mass1,invMass1,Mass01,Mass2,invMass2,Mass02,Mass12]=MassLGTriangle();
%
Uxex=zeros(NdofU,1); Uyex=Uxex; 
%
RHS_Uxex=Uxex; RHS_Uyex=RHS_Uxex;
RHS_Ux=RHS_Uxex; RHS_Uy=RHS_Uyex;
%
if (ordre==2)
  ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
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
  UxexLOC=awp.*(1-cosxy(:,1)).*sinxy(:,2);
  UyexLOC=awp.*(cosxy(:,2)-1).*sinxy(:,1);
  if (ordre==1)
    psiLG=lambda; UGLO=IGLO;
  end
  %
  if (ordre==2)
    AGLO=TriEdg(t,:);
    lambda2=lambda.*lambda;
    psiS=2*lambda2-lambda;
    psiM(:,ILOC)=4*lambda(:,JLOC).*lambda(:,KLOC);
    psiLG=[psiS,psiM]; UGLO=[IGLO,AGLO+Nbpt]; 
  endif
  RHS_Ux(UGLO)+=sum(UxexLOC.*psiLG,1)';
  RHS_Uy(UGLO)+=sum(UyexLOC.*psiLG,1)';
end
%
%%%%%%%%%%%%%%
% Vitesse
%%%%%%%%%%%%%%
Uxex=Mu\RHS_Ux; Uyex=Mu\RHS_Uy;
