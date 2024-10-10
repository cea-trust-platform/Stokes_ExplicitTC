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
% ErrorStokesSinusSV.m:
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
% SYNOPSIS [Eu0,Eu1,Ep0]=ErrorStokesSinusSV(Mu,Ku,Uxh,Uyh,UxP,UyP,UxI,UxI,Mp,Ph,PeP,PeI,moyPeP)
% GLOBAL - Nbpt : nombre de sommets
%        - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - Nbtri : nombre de triangles
%        - Aires : aire des triangles
%        - ordre : 1 ou 2
%
% OUTPUT - Eu0 : erreur L2 vitesse
%        - Eu1 : erreur H1 vitesse
%        - Ep0 : erreur L2 pression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Eu0,Eu1,Ep0]=ErrorStokesSinusSV(Mu,Ku,Uxh,Uyh,UxP,UyP,UxI,UyI,Mp,Ph,PeP,PeI,moyPexP)
global Nbpt CoorNeu
global Nbtri Aires NumTri 
global Nbedg TriEdg EdgNorm EdgTri
global ordre
global nu npi npi2 bp
%
nu2=nu*nu;
normPex2=0.25*bp^2;
normGUex2=2*npi2;
normUP2=normGUex2+normPex2/nu2;
nu2normUP2=nu2*normGUex2+normPex2;
%
nErr=3;
Eu0=zeros(1,nErr); Eu1=zeros(1,nErr); Ep0=zeros(1,nErr);
%
dUxhP=UxP-Uxh; dUyhP=UyP-Uyh; dUxhI=UxI-Uxh; dUyhI=UyI-Uyh;
Eu0(1)=dUxhP'*Mu*dUxhP+dUyhP'*Mu*dUyhP; Eu0(2)=dUxhI'*Mu*dUxhI+dUyhI'*Mu*dUyhI;
Eu1(1)=dUxhP'*Ku*dUxhP+dUyhP'*Ku*dUyhP; Eu1(2)=dUxhI'*Ku*dUxhI+dUyhI'*Ku*dUyhI;
dPhP=PeP-Ph; dPhI=PeI-Ph;
Ep0(1)=dPhP'*Mp*dPhP; Ep0(2)=dPhI'*Mp*dPhI;
%
npi2nu=npi2*nu;
cx=npi*bp-2*npi2nu;
cy=npi*bp+2*npi2nu;
%
ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
for t=1:Nbtri
  IGLO=NumTri(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  % volume of the element
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp; coeff=1/(2*aire);
  sinxy=sin(npi*xyp); cosxy=cos(npi*xyp);
  UxLOC = (1-cosxy(:,1)).*sinxy(:,2);
  GxUxLOC=npi*sinxy(:,1).*sinxy(:,2); GyUxLOC=npi*(1-cosxy(:,1)).*cosxy(:,2);
  UyLOC = (cosxy(:,2)-1).*sinxy(:,1);
  GxUyLOC=npi*(1-cosxy(:,2)).*cosxy(:,1); GyUyLOC=-npi*sinxy(:,1).*sinxy(:,1);
  PexLOC=bp*sinxy(:,1).*sinxy(:,2);
  % vecteurs face-normale
  AGLO=TriEdg(t,:); EdgNormT=EdgNorm(AGLO,:);
  for iloc=1:2
    if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
    end
  end
  EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
  GxLambda=-coeff*(EdgNormT(:,1)*ones(1,7))';
  GyLambda=-coeff*(EdgNormT(:,2)*ones(1,7))';
  % EF P1-P0 (ordre=1)
  if (ordre==1)
    phiP=1;
    phiLG=lambda; UGLO=IGLO;
    GxPhiLG=GxLambda;
    GyPhiLG=GyLambda;
    PhT=Ph(t,1);
  end
  % EF P2-P1dg (ordre=2)
  if (ordre==2)
    phiP=lambda;
    lambda2=lambda.*lambda;
    phiS=2*lambda2-lambda;
    phiM(:,ILOC)=4*lambda(:,JLOC).*lambda(:,KLOC);
    phiLG=[phiS,phiM]; UGLO=[IGLO,AGLO+Nbpt];
    qL=4*lambda; qLm1=qL-1;
    GxPhiLG=[qLm1.*GxLambda,qL(:,JLOC).*GxLambda(:,KLOC)+qL(:,KLOC).*GxLambda(:,JLOC)];
    GyPhiLG=[qLm1.*GyLambda,qL(:,JLOC).*GyLambda(:,KLOC)+qL(:,KLOC).*GyLambda(:,JLOC)];
    fin=3*t;
    PhT=Ph(fin-2:fin,1);
  end
  %
  UxhLOC=phiLG*Uxh(UGLO,1); UyhLOC=phiLG*Uyh(UGLO,1);
  GxUxhLOC=GxPhiLG*Uxh(UGLO,1); GyUxhLOC=GyPhiLG*Uxh(UGLO,1);
  GxUyhLOC=GxPhiLG*Uyh(UGLO,1); GyUyhLOC=GyPhiLG*Uyh(UGLO,1);
  PhLOC=phiP*PhT;
  %
  Eu0(3)+=awp*((UxLOC-UxhLOC).^2+(UyLOC-UyhLOC).^2);
  Eu1(3)+=awp*((GxUxLOC-GxUxhLOC).^2+(GyUxLOC-GyUxhLOC).^2+(GxUyLOC-GxUyhLOC).^2+(GyUyLOC-GyUyhLOC).^2);
  Ep0(3)+=awp*((PexLOC-PhLOC).^2);
  %
end
%
Eu0=sqrt(Eu0/normUP2);
Eu1=sqrt(Eu1/normUP2);
Ep0=sqrt(Ep0/nu2normUP2);