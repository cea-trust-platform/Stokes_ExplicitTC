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
% SolutionStokesSinusNC.m:
% Solution du probleme de Stokes "Sinus"
%   avec les elements finis non conformes CR ou FS 
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= npi*[ sin(npi*y)*{(b-2*nu*npi)cos(npi*x)+nu*npi} ; 
%                                sin(npi*x)*{(b+2*nu*npi)cos(npi*y)-nu*npi)}]
%    div U = 0;
%
% U(x,y)=[(1-cos(npi*x))*sin(npi*y) ; 
%         (cos(npi*y)-1)*sin(npi*x)]; 
% p(x,y)=sin(npi*x)*sin(npi*y);
%
%
% SYNOPSIS [Uxex,Uyx,Pex,Uex2h,GUex2h,Pex2h]=SolutionStokesSinusNC(Mu,Ku,Mp)
%
% GLOBAL - Nbtri            : nb de triangles
%        - CoorNeu(Nbpt,2)  : coordonnees (x, y) des sommets
%        - CoorMil(Nbedg,2) : coordonnees (x, y) des milieux d'aretes
%        - CoorBary(Nbtri,2): coordonnees (x, y) des barycentres des triangles
%        - NumTri(Nbtri,3)  : liste de triangles (3 numeros de sommets)
%        - Nbedg            : nb d'aretes
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
% OUTPUT - [Uxex,Uyex] : la solution exacte au points de discretisation
%        - Uex2h,GUex2h : les normes de la solution projetee 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxex,Uyex,Uex2h,GUex2h]=SolutionStokesSinusNC(Mu,Ku)
global Nbpt CoorNeu RefNeu
global Nbtri CoorBary NumTri Aires
global Nbedg CoorMil TriEdg RefEdg NumEdg
global ordre npi
%
NdofU=Nbedg;
UxexF=zeros(Nbedg,1); UyexF=UxexF;
if (ordre==2)
  NP2=Nbedg+Nbpt;
  NdofU=NP2+Nbtri;
end
Uxex=zeros(NdofU,1); Uyex=Uxex;
%
for e=1:Nbedg
  IGLO=NumEdg(e,:);
  CoorNeuE=CoorNeu(IGLO,:);
  [xyp,wp,lambda,np]=IntEdg_Boo5(CoorNeuE);
  npixyp=npi*xyp;
  cosxy=cos(npixyp); sinxy=sin(npixyp);
  UxexF(e,1)=wp*((1-cosxy(:,1)).*sinxy(:,2));
  UyexF(e,1)=wp*((cosxy(:,2)-1).*sinxy(:,1));
endfor
%
if (ordre==1)
  Uxex=UxexF;
  Uyex=UyexF;
end
%
if (ordre==2)
  npiXY=npi*CoorNeu;
  cosXY=cos(npiXY); sinXY=sin(npiXY);
  DofU_S=find(RefNeu==0);
  Uxex(DofU_S,1)=(1-cosXY(DofU_S,1)).*sinXY(DofU_S,2);
  Uyex(DofU_S,1)=(cosXY(DofU_S,2)-1).*sinXY(DofU_S,1);
  %
  for e=1:Nbedg
    IGLO=NumEdg(e,:);
    Uxex(Nbpt+e,1)=1.5*UxexF(e,1)-0.25*sum(Uxex(IGLO,1));
    Uyex(Nbpt+e,1)=1.5*UyexF(e,1)-0.25*sum(Uyex(IGLO,1));
  end
  %
end
%
GUex2h=Uxex'*Ku*Uxex+Uyex'*Ku*Uyex;
Uex2h=Uxex'*Mu*Uxex+Uyex'*Mu*Uyex;