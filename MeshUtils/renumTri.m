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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot CEA
%
% renumTri.m:
% routine de renumerotation des triangles pour DGFEM
% 
% SYNOPSIS  [newnum,NumTriNew,RefTriNew,NumTri2New,RefTri2New,CoorMilNew,AiresNew,EdgTriNew]=...
%             renumTri(NumTri,NumTri2,RefTri,RefTri2,TriEdg,CoorMil,Aires,EdgTri,RefEdg)
%          
% INPUT  - Donnees du maillage, numerotation initiale
%
% OUTPUT - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - RefTri(Nbtri,1) : Reference de chaque triangle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newnum,NumTriNew,RefTriNew,NumTri2New,RefTri2New,TriEdgNew,CoorMilNew,AiresNew,EdgTriNew]=...
         renumTri(NumTri,NumTri2,RefTri,RefTri2,TriEdg,CoorMil,Aires,EdgTri,RefEdg)
  Nbtri=size(Numtri,1);
  EdgI=find(RefEdg==0);
  EdgTrI=EdgTri(EdgI);
  NbedgI=size(EdgTriI,1);
  Pro=speye(Nbtri,Nbtri);
  for e=1:NbedgI
    T1=EdgTrI(e,1);
    T2=EdgTrI(e,2);
    Pro(T1,T2)=1;
  end
  newnum=symrcm(Pro);
  NumTriNew=NumTri(newnum,:);
  RefTriNew=RefTri(newnum,:);
  TriEdgNew=TriEdg(newnum,:);
  CoorMilNew=CoorMil(newnum,:);
  AiresNew=Aires(newnum,:);
  Nbedg=size(EdgTri,1);
  EdgTriNew=EdgTri;
  for e=1:Nbedg
    T1=EdgTri(e,1);
    EdgTriNew(e,1)=newnum(T1);
    T2=EdgTri(e,2);
    if (T2>0)
      EdgTriNew(e,2)=newnum(T2);
    end   
  end
  NumTri2New=zeros(4*Nbtri,3);
  RefTri2New=zeros(4*Nbtri,1);
  for tt=1:Nbtri
    fin=4*tt;
    deb=fin-3;
    tnew=newnum(tt);
    finnew=4*tnew;
    debnew=fin-3;
    NumTri2New(debnew:finnew,:)=NumTri2(deb:fin,:);
    RefTri2New(debnew:finnew,:)=RefTri(deb:fin,:);
  end