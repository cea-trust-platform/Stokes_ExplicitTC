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
% writepolco.m:
% routine d'ecriture de fichiers de maillages triangulaires 2D au format .edg 
% 
% SYNOPSIS writepolco(rho2,rho,theta,xy0tri,xy0edg,xy0som)
%          
% INPUT  : - rho2   : distance des points a xy0 au carre
%          - rho    : distance des points a xy0
%          - theta  : angle polaire dans le repere (xy0,X,Y)
%          - xy0tri : triangles ayant xy0 comme sommet
%          - xy0edg : aretes ayant xy0 comme sommet
%          - xy0som : numero du sommet sur lequel se trouve xy0
%          - filename : le nom d'un fichier de maillage SANS SUFFIXE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writepolco(rho2,rho,theta,xy0tri,xy0edg,xy0som,filename)
  % OK pour maillage d'un carre [0,1]x[0,1], 
  % OK pour xy0 en (0.5,0.5) correspondant a un sommet
  polfile=strcat(filename,".pol");
  fid2=fopen(polfile,'w');
  Ndof=size(rho2,1);
  Ntri=size(xy0tri,1);
  Nedg=size(xy0edg,1);
  %  
  fprintf(fid2,'$Nxy0som\n');
  fprintf(fid2,'%i\n',xy0som);
  fprintf(fid2,'$EndNxy0som\n');
  %
  fprintf(fid2,'$Ntri\n');
  fprintf(fid2,'%i\n',Ntri);
  for t=1:Ntri
    fprintf(fid2,'%i %i\n',xy0tri(t,:));
  end
  fprintf(fid2,'$EndNtri\n');
  fprintf(fid2,'$Nedg\n');
  fprintf(fid2,'%i\n',Nedg);
  for e=1:Nedg
    fprintf(fid2,'%i %i\n',xy0edg(e,:));
  end
  fprintf(fid2,'$EndNedg\n');
  fprintf(fid2,'$Ndof\n');
  fprintf(fid2,'%i\n',Ndof);
  %
  for n=1:Ndof
    fprintf(fid2,'%i %2.18e %2.18e %2.18e\n',n,rho2(n),rho(n),theta(n));
  end
  fprintf(fid2,'$EndNdof\n');
  fclose(fid2);
  
  
  