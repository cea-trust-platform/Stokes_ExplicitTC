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
% buildSVmesh.m:
% Raffinement d'un maillage P1 avec 3 sous-triangles par triangles en utilisant les barycentres
% Utile pour la discretisation Scott-Vogelius, P2-P1disc
%
% SYNOPSIS [CoorNeuSV,CoorNeu2SV,CoorBarySV,RefNeuSV,RefNeu2SV,NumTriSV,NumTri2SV,RefTriSV,RefTri2SV,NumEdgSV,CoorMilSV,...
%           RefEdgSV,TriEdgSV,EdgTriSV,SomOppSV,Lga2SV,EdgNormSV,AiresSV]=...
%           buildSVmesh(CoorNeu,CoorBary,RefNeu,NumTri,RefTri,NumEdg,CoorMil,RefEdg,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires)
%          
% INPUT  - CoorNeu(Nbpt,2)   : coordonnees (x, y) des sommets
%        - CoorBary(Nbtri,2) : coordonnees (x, y) des barycentres des elements
%        - RefNeu(Nbpt,1)    : reference des sommets
%        - NumTri(Nbtri,3)   : liste de triangles (3 numeros de sommets)
%        - RefTri(Nbtri,1)   : reference de chaque triangle
%        - NumEdg(NbEdg,2)   : liste d'aretes (2 numeros de sommets)
%        - CoorMil(NbEdg,2)  : coordonnees (x, y) des milieux d'aretes
%        - RefEdg(NbEdg,1)   : reference de chaque arete 
%        - TriEdg(Nbtri,3)   : Pour chaque triangle l, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%        - EdgTri(NbEdg,2)   : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere 
%        - SomOpp(NbEdg,2)   : Numero du sommet oppose a l'arete dans chaque triangle
%                                  SomOpp(a,2) = 0 si a est sur la frontiere 
%        - Lga2(NbEdg,1)     : longueurs des aretes au carre
%        - EdgNorm(NbEdg,2)  : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1)    : aires des triangles
%
% OUTPUT - CoorNeuSV(NbptSV,2)   : coordonnees (x, y) des sommets
%        - CoorNeu2SV(NbptSV2,2) : coordonnees (x, y) des sommets suivis des milieux d'aretes
%        - CoorBarySV(NbtriSV,2) : coordonnees (x, y) des barycentres des elements
%        - RefNeuSV(NbptSV,1)    : reference des sommets
%        - RefNeu2SV(NbptSV2,1)  : reference des sommets suivis des milieux d'aretes
%        - NumTriSV(NbtriSV,3)   : liste de triangles (3 numeros de sommets)
%        - NumTri2SV(NbtriSV2,3) : liste de triangles du maillage raffine (3 numeros de noeuds)
%        - RefTriSV(NbtriSV,1)   : reference de chaque triangle
%        - RefTri2SV(NbtriSV2,1) : reference des triangles du maillage raffine
%        - NumEdgSV(NbEdgSV,2)   : liste d'aretes (2 numeros de sommets)
%        - CoorMilSV(NbEdgSV,2)  : coordonnees (x, y) des milieux d'aretes
%		     - RefEdgSV(NbEdgSV,1)   : reference de chaque arete 
%		     - TriEdgSV(NbtriSV,3)   : Pour chaque triangle l, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTriSV(NbEdgSV,2)   : Pour chaque arete, EdgTriSV(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTriSV(a,2) = 0 si a est sur la frontiere 
%		     - SomOppSV(NbEdgSV,2)   : Numero du sommet oppose a l'arete dans chaque triangle
%                                  SomOppSV(a,2) = 0 si a est sur la frontiere 
%        - Lga2SV(NbEdgSV,1)     : longueurs des aretes au carre
%        - EdgNormSV(NbEdgSV,2)  : vecteurs face-normale, orientes tri1->tri2
%        - AiresSV(NbtriSV,1)    : aires des triangles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoorNeuSV,CoorNeu2SV,CoorBarySV,RefNeuSV,RefNeu2SV,NumTriSV,NumTri2SV,RefTriSV,RefTri2SV,NumEdgSV,CoorMilSV,...
           RefEdgSV,TriEdgSV,EdgTriSV,SomOppSV,Lga2SV,EdgNormSV,AiresSV]=...
           buildSVmesh(CoorNeu,CoorBary,RefNeu,NumTri,RefTri,NumEdg,CoorMil,RefEdg,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires)
%
Nbpt=size(CoorNeu,1);
Nbtri=size(NumTri,1);
NbEdg=size(EdgTri,1);
%
CoorNeuSV  = [CoorNeu;CoorBary];
CoorBarySV = zeros(3*Nbtri,2);
RefNeuSV   = [RefNeu;zeros(Nbtri,1)];
NumTriSV   = zeros(3*Nbtri,3);
RefTriSV   = zeros(3*Nbtri,1);
NumEdgSV   = [NumEdg ;zeros(3*Nbtri,2)];
CoorMilSV  = [CoorMil;zeros(3*Nbtri,2)];
RefEdgSV   = [RefEdg ;zeros(3*Nbtri,1)];
TriEdgSV   = zeros(3*Nbtri,3);;
EdgTriSV   = zeros(NbEdg+3*Nbtri,2);
SomOppSV   = zeros(NbEdg+3*Nbtri,2);
Lga2SV     = [Lga2;zeros(3*Nbtri,1)];
EdgNormSV  = [EdgNorm;zeros(3*Nbtri,2)];
AiresSV    = zeros(3*Nbtri,1);
% A incrementer
NbEdgSV= NbEdg;
NbptSV=Nbpt;
NbTriSV=0;
JK=[2,3,1;
    3,1,2];
for t=1:Nbtri
  aire3=Aires(t)/3;
  % Nouveau noeud
  NbptSV=NbptSV+1;
  for iloc=1:3
    jloc=JK(1,iloc); 
    kloc=JK(2,iloc);
    % Ancienne arete
    SjSk=TriEdg(t,iloc);
    NbTriT=NbTriSV+iloc;
    for tt=1:2
        if (EdgTri(SjSk,tt)==t)
           EdgTriSV(SjSk,tt)=NbTriT;
           SomOppSV(SjSk,tt)=NbptSV;
        end
    end
    % Nouveau triangle
    NumTriSV(NbTriT,:)=[NumTri(t,jloc),NumTri(t,kloc),NbptSV];
    TriEdgSV(NbTriT,:)=[ NbEdgSV+kloc , NbEdgSV+jloc ,SjSk ];
    RefTriSV(NbTriT,1)=RefTri(t,1);
    AiresSV (NbTriT,1)=aire3;
    % Nouvelle arete
    NbEdgT=NbEdgSV+iloc;
    EdgTriSV (NbEdgT,:)=[NbTriSV+jloc,NbTriSV+kloc];
    Si=NumTri(t,iloc);
    CoorMilSV(NbEdgT,:)=0.5*(CoorBary(t,:)+CoorNeu(Si,:));
    SomOppSV (NbEdgT,:)=[NumTri(t,kloc),NumTri(t,jloc)];
    NumEdgSV (NbEdgT,:)=[Si,NbptSV];
  end  
  % Anciennes faces normales
  ALOC=TriEdg(t,:);
  EdgNormT=EdgNorm(ALOC,:);
  for i=1:2
    a=TriEdg(t,i);  
    if (EdgTri(a,1)~=t)
       EdgNormT(i,:)=-EdgNormT(i,:);
    end       
  end
  EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
  % Nouvelles faces normales
  E1=NbEdgSV+1; E2=NbEdgSV+2; E3=NbEdgSV+3;
  S1=NumTri(t,1); S3=NumTri(t,3);
  EdgNormSV(E1,:)=[CoorBary(t,2)-CoorNeu(S1,2),CoorNeu(S1,1)-CoorBary(t,1)];
  S3M1=CoorMilSV(E1,:)-CoorNeu(S3,:);
  % EdgNormE1 doit sortir de T2
  if (dot(EdgNormSV(E1,:),S3M1)<0)
    EdgNormSV(E1,:)=-EdgNormSV(E1,:);
  end
  Lga2SV(E1)=dot(EdgNormSV(E1,:),EdgNormSV(E1,:));
  % EdgNormE2 doit sortir de T3
  EdgNormSV(E2,:)= EdgNormSV(E1,:)-EdgNormT(3,:); Lga2SV(E2,1)=dot(EdgNormSV(E2,:),EdgNormSV(E2,:));
  % EdgNormE3 doit sortir de T1
  EdgNormSV(E3,:)= EdgNormSV(E2,:)-EdgNormT(1,:); Lga2SV(E3,1)=dot(EdgNormSV(E3,:),EdgNormSV(E3,:));
  % On incremente
  NbTriSV=NbTriT;
  NbEdgSV=NbEdgT;
end

% VERIFICATIONS
for a=1:NbEdgSV
  EdgNorm_a=EdgNormSV(a,:);
  T1=EdgTriSV(a,1);
  S1=SomOppSV(a,1);
  S1M1=CoorMilSV(a,:)-CoorNeuSV(S1,:);
  if (dot(EdgNorm_a,S1M1)<0)
     dot(EdgNorm_a,S1M1)
     printf('Arete SV %i dans le mauvais sens\n',a);
  end   
end
for t=1:NbTriSV
  X=zeros(3,3);
  for iloc=1:3
      i=NumTriSV(t,iloc);
      X(iloc,1:2)=CoorNeuSV(i,:);
  end
  % sens trigo ?
  S1S2=[X(2,:)-X(1,:)];
  S1S3=[X(3,:)-X(1,:)];
  z=cross(S1S2,S1S3);
  sens=1;
  if (z(3)<0)
     printf('Triangle SV %i sens non trigo\n',t);
     sens=-1;
     % on echange les deux derniers elements
     tmpT=NumTriSV(t,2);       
     NumTriSV(t,2)=NumTriSV(t,3);
     NumTriSV(t,3)=tmpT;
     %
     tmpE=TriEdgSV(t,2);
     NumTriSV(t,2)=TriEdgSV(t,3);
     TriEdgSV(t,3)=tmpE;
  end
end
[CoorNeu2SV,RefNeu2SV,NumTri2SV,RefTri2SV]=buildP2mesh(CoorNeuSV,RefNeuSV,NumTriSV,RefTriSV,NumEdgSV,CoorMilSV,RefEdgSV,TriEdgSV);