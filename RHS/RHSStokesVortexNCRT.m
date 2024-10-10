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
% RHSStokesVortexNCRT.m:
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
% SYNOPSIS [RHS_Ux,RHS_Uy]=RHSStokesVortexNCRT()
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
%        - RHS_p  : second membre pression (relevement CLDNH)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHS_Ux,RHS_Uy]=RHSStokesVortexNCRT( )
global Nbpt CoorNeu RefNeu
global Nbtri Aires NumTri CoorBary RefTri SomOpp
global Nbedg TriEdg RefEdg EdgTri EdgNorm NumEdg LgEdg2 CoorMil
global EdgC xy0som xy0tri
global nu alpha beta ordre
%
if (ordre==1)
  basis='CRRT';
end
%
if (ordre==2)
  basis='FSRT';
end
%
NdofU=Nbedg;
if (ordre==2)
  NP2=Nbedg+Nbpt;
  NdofU=NP2+Nbtri;
end
RHS_Ux=zeros(NdofU,1); RHS_Uy=zeros(NdofU,1);
UxexNC=zeros(NdofU,1); UyexNC=zeros(NdofU,1);
%
ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
cnu=-nu*(alpha^2-1);
%
% Triangles qui touchent le centre
[RHS_Ux,RHS_Uy]=RHSStokesVortexPolTri(nu,alpha,beta,basis,NdofU);
%
% Triangles qui ne touchent pas le centre
TriE=find(RefTri==0); NbtriE=size(TriE,1);
% EF de CR P1-NC (ordre=1)
if (ordre==1)
  for tri=1:NbtriE
    t=TriE(tri,1);
  %for t=1:Nbtri
    IGLO=NumTri(t,:);
    CoorNeuT=CoorNeu(IGLO,:);
    % Integration points
    [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
    awp=0.5*wp';
    X0=xyp(:,1)-0.5; Y0=xyp(:,2)-0.5;
    theta=atan2(Y0,X0);
    rho=sqrt(X0.^2+Y0.^2);
    cth=cos(theta); sth=sin(theta);
    rhoA=rho.^(alpha-2); rhoB=rho.^(beta-1);
    FxLOC = awp.*(-cnu*rhoA.*sth+beta*rhoB.*cth);
    FyLOC = awp.*( cnu*rhoA.*cth+beta*rhoB.*sth);
    % (d|T|)^{-1}\int_T(x-OSi)\cdot\Fvec
    RHS_tmp=(ones(3,1)*xyp(:,1)'-CoorNeuT(:,1)*ones(1,np))*FxLOC+...
            (ones(3,1)*xyp(:,2)'-CoorNeuT(:,2)*ones(1,np))*FyLOC;
    %
    % vecteurs face-normale
    AGLO=TriEdg(t,:); EdgNormT=EdgNorm(AGLO,:);
    for iloc=1:2
     if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
     end
    end
    EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
    RHS_Ux(AGLO)+=EdgNormT(:,1).*RHS_tmp;
    RHS_Uy(AGLO)+=EdgNormT(:,2).*RHS_tmp;
    %
  end
end
%
% EF de FS (ordre=2) 
if (ordre==2)
  ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
  phiM=zeros(7,3);
  for tri=1:NbtriE
    t=TriE(tri,1);
    IGLO=NumTri(t,:);
    CoorNeuT =CoorNeu(IGLO,:);
    % vecteurs face-normale
    AGLO=TriEdg(t,:); EdgNormT=EdgNorm(AGLO,:);
    for iloc=1:2
     if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
     end
    end
    EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
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
    %
    [PiRTxx,PiRTxy,PiRTyx,PiRTyy]=PiRT1(aire,CoorNeuT,EdgNormT,lambda);
    %
    UGLO=[IGLO,AGLO+Nbpt,NP2+t];
    RHS_Ux(UGLO)+=(FxLOC'*PiRTxx+FyLOC'*PiRTxy)';
    RHS_Uy(UGLO)+=(FxLOC'*PiRTyx+FyLOC'*PiRTyy)';
    % 
  end
end
%