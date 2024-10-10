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
% writeSVmesh.m:
% routine d'ecriture de fichiers de maillages triangulaires 2D obtenus par raffinement barycentrique
% 
% SYNOPSIS writeSVmesh(NumEdgB,RefEdgB,CoorNeuSV,CoorBarySV,RefNeuSV,NumTriSV,RefTriSV,NumEdgSV,CoorMilSV,RefEdgSV,TriEdgSV,EdgTriSV,...
%                     SomOppSV,Lga2SV,EdgNormSV,AiresSV,filename)
%          
% INPUT  - NumEdgB(NbEdgB,1)   : liste des aretes du bord (2 numeros de sommets)
%        - RefEdgB(NbEdgB,1)   : reference des aretes du bord
%        - CoorNeuSV(NbptSV,2) : coordonnees (x, y) des sommets maillage SV
%        - CoorNeu2SV(Nbpt2,2) : coordonnees (x, y) des sommets suivis des milieux d'aretes
%        - CoorBarySV(Nbtri,2) : coordonnees (x, y) des barycentres des elements
%        - RefNeuSV(Nbpt,1)    : reference des sommets
%        - RefNeu2SV(Nbpt2,1)  : reference des sommets suivis des milieux d'aretes
%        - NumTriSV(Nbtri,3)   : liste de triangles (3 numeros de sommets)
%        - NumTri2SV(Nbtri2,3) : liste de triangles du maillage raffine (3 numeros de noeuds)
%        - RefTriSV(Nbtri,1)   : reference de chaque triangle
%        - RefTri2SV(Nbtri2,1) : reference des triangles du maillage raffine
%        - NumEdgSV(NbEdgSV,2) : liste d'aretes (2 numeros de sommets)
%        - CoorMilSV(NbEdgSV,2): coordonnees (x, y) des milieux d'aretes
%		     - RefEdgSV(NbEdgSV,1) : reference de chaque arete 
%		     - TriEdgSV(NbtriSV,3) : Pour chaque triangle l, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTriSV(NbEdgSV,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTriSV(a,2) = 0 si a est sur la frontiere 
%		     - SomOppSV(NbEdgSV,2) : Numero du sommet oppose a l'arete dans chaque triangle
%                                  SommetOppose(a,2) = 0 si a est sur la frontiere 
%        - Lga2SV(NbEdgSV,1)   : longueurs des aretes au carre
%        - EdgNormSV(NbEdgSV,2): vecteurs face-normale, orientes tri1->tri2
%        - AiresSV(NbtriSV,1)  : aires des triangles
%        - filename : le nom d'un fichier de maillage de sauvegarde
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SVfilename=writeSVmesh(NumEdgB,RefEdgB,CoorNeuSV,CoorNeu2SV,CoorBarySV,RefNeuSV,NumTriSV,NumTri2SV,RefTriSV,RefTri2SV,...
                     RefNeu2SV,NumEdgSV,CoorMilSV,RefEdgSV,TriEdgSV,EdgTriSV,...
                     SomOppSV,Lga2SV,EdgNormSV,AiresSV,filename)
% Enregistrement des donnees maillage SV.
SVfilename=strcat(filename,"_SV");
meshfile=strcat(SVfilename,".msh");
fid=fopen(meshfile,'w');
%
fprintf(fid,'$Nodes\n');
    
NbptSV=size(CoorNeuSV,1);
fprintf(fid,'%i\n',NbptSV);
z=0;
for i=1:NbptSV
    fprintf(fid,'%i %2.18e %2.18e %i\n',i,CoorNeuSV(i,:),z);
end
fprintf(fid,'$EndNodes\n');
fprintf(fid,'$Elements\n');
NbtriSV=size(NumTriSV,1);
NbEdgB=size(NumEdgB,1);
NbelemSV=NbtriSV+NbEdgB;
fprintf(fid,'%i\n',NbelemSV);
dim1=1;
dim2=2;
for a=1:NbEdgB
  fprintf(fid,'%i %i %i %i %i %i %i\n',a,dim1,dim2,RefEdgB(a,1),RefEdgB(a,1),NumEdgB(a,:));
end
for t=1:NbtriSV
  fprintf(fid,'%i %i %i %i %i %i %i %i\n',t+NbEdgB,dim2,dim2,RefTriSV(t,1),dim1,NumTriSV(t,:));
end
fprintf(fid,'$EndElements\n');
fclose(fid);
% Enregistrement des donnees des aretes.
writeedges(CoorBarySV,NumEdgSV,CoorMilSV,RefEdgSV,TriEdgSV,EdgTriSV,SomOppSV,Lga2SV,EdgNormSV,AiresSV,SVfilename);
% Enregistrement des donnees P2.
writeP2mesh(CoorNeu2SV,RefNeu2SV,NumTri2SV,RefTri2SV,SVfilename);