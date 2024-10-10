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
% readpolco.m:
% routine de lecture des fichiers de maillages triangulaires 2D au format .pol 
% 
% SYNOPSIS [rho2,rho,theta,xy0tri,xy0edg,xy0som]=readpolco(filename)
%          
% INPUT  - filename : le nom d'un fichier au format pol sans son suffixe .pol
%
% OUTPUT : - rho2   : distance des points a xy0 au carre
%          - rho    : distance des points a xy0
%          - theta  : angle polaire dans le repere (xy0,X,Y)
%          - xy0tri : triangles ayant xy0 comme sommet
%          - xy0edg : arete sur laquelle se trouve xy0
%          - xy0som : sommet sur lequel se trouve xy0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho2,rho,theta,xy0tri,xy0edg,xy0som]=readpolco(filename)
polcofile=strcat(filename,".pol");
fid2=fopen(polcofile,'r');
if fid2 <=0,
  msg=['Le fichier de maillage : ' polcofile ' n''a pas ete trouve'];
  error(msg);
end
%% Sommet
while ~strcmp(fgetl(fid2),'$Nxy0som'), end
xy0som = str2num(fgetl(fid2));
%
%% Triangles
%
while ~strcmp(fgetl(fid2),'$Ntri'), end
Ntri= str2num(fgetl(fid2));
xy0tri=zeros(Ntri,2);
% boucle sur les Tri
for t=1:Ntri
  xy0tri(t,:)= str2num(fgetl(fid2));
end
while ~strcmp(fgetl(fid2),'$Nedg'), end
%
%% Edges
%
Nedg= str2num(fgetl(fid2));
xy0edg=zeros(Nedg,2);
% boucle sur les Edg
for t=1:Nedg
  xy0edg(t,:)= str2num(fgetl(fid2));
end
while ~strcmp(fgetl(fid2),'$Ndof'), end
Ndof= str2num(fgetl(fid2));
rho2=zeros(Ndof,1); rho=rho2; theta=rho2;
% boucle sur les degres de liberte
for n=1:Ndof
  tmp= str2num(fgetl(fid2));
  rho2(n) = tmp(2);
  rho(n)  = tmp(3);
  theta(n)= tmp(4);
end
