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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot CEA
%
% buildmeshfiles.m:
% routine de lecture de fichiers de maillages triangulaires 2D au format .msh
% 
% SYNOPSIS [CoorNeu,RefNeu,NumTri,RefTri,NumEdgB,RefEdgB]=buildmeshfiles(filename)
%          
% INPUT  - filename : le nom d'un fichier de maillage au format msh
%                   SANS SON SUFFIXE .msh
%
% OUTPUT - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - RefTri(Nbtri,1) : Reference de chaque triangle
%        - NumEdgB(NbEdgB,2) : Numero des 2 noeuds de chaque arete
%		     - RefEdgB(NbEdgB,1) : Reference de chaque arete 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buildmeshfiles(filename)
  % Lecture du maillage initial
  [CoorNeu,RefNeu,NumTri,RefTri,NumEdgB,RefEdgB]=readmesh(filename);
  % Calcul des normales aux faces 
  [CoorBary,NumEdg,CoorMil,RefEdg,TriEdg,...
          EdgTri,SomOpp,Lga2,EdgNorm,Aires]=buildedges(CoorNeu,RefNeu,NumTri,NumEdgB,RefEdgB);
  % Stockage des informations relatives aux aretes
  writeedges(CoorBary,NumEdg,CoorMil,RefEdg,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires,filename);
  % Construction maillage P2
  [CoorNeu2,RefNeu2,NumTri2,RefTri2]=buildP2mesh(CoorNeu,RefNeu,NumTri,RefTri,NumEdg,CoorMil,RefEdg,TriEdg);
  % Stockage des informations relatives au maillage P2
  writeP2mesh(CoorNeu2,RefNeu2,NumTri2,RefTri2,filename);
  % construction maillage de Scott-Vogelius
  [CoorNeuSV,CoorNeu2SV,CoorBarySV,RefNeuSV,RefNeu2SV,NumTriSV,NumTri2SV,RefTriSV,RefTri2SV,NumEdgSV,CoorMilSV,...
           RefEdgSV,TriEdgSV,EdgTriSV,SomOppSV,Lga2SV,EdgNormSV,AiresSV]=buildSVmesh(CoorNeu,CoorBary,RefNeu,NumTri,RefTri,NumEdg,CoorMil,...
           RefEdg,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires);
  % Stockage, ecriture maillage Scott-Vogelius
  SVfilename=writeSVmesh(NumEdgB,RefEdgB,CoorNeuSV,CoorNeu2SV,CoorBarySV,RefNeuSV,NumTriSV,NumTri2SV,RefTriSV,RefTri2SV,...
                     RefNeu2SV,NumEdgSV,CoorMilSV,RefEdgSV,TriEdgSV,EdgTriSV,...
                     SomOppSV,Lga2SV,EdgNormSV,AiresSV,filename);