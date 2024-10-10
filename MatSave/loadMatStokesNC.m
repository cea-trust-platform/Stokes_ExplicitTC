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
% loadMatStokesNC.m:
%
% Matrices de raideur et matrice de masse
% Formulation DG SIP, fonctions de base monomiales
%
% SYNOPSIS [Ku,Mu,Mp,invMp,Bx,By]=loadMatStokesNC(mi)
%
%
% OUTPUT - Ku : matrice de raideur vitesse
%        - Mu : matrice de masse vitesse
%        - Mp : matrice de masse pression
%        - invMp : inverse de la matrice de masse pression
%        - Bx : matrice b(vx,q)
%        - By : matrice b(vy,q)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ku,Mu,Mp,invMp,Bx,By]=loadMatStokesNC(filename)

% VARIABLES GLOBALES
global ordre
global Nbpt
global Nbtri Aires
global Nbedg RefEdg EdgTri
%
NP1=3*Nbtri; NP2=6*Nbtri;
%
if (ordre==1)
  Nu=Nbedg; Np=Nbtri;
end
%
if (ordre==2)
  Nu=Nbpt+Nbedg+Nbtri; Np=3*Nbtri;
end
%
oU=ordre; oP=ordre-1;
%  
Kufile=sprintf('KuNC%i_%s',oU,filename);
Ku=loadMat(Kufile,Nu,Nu);
%
Mufile=sprintf('MuNC%i_%s',oU,filename);
Mu=loadMat(Mufile,Nu,Nu);
%
if (ordre==1)
  Mp=diag(sparse(Aires));
  invMp=diag(sparse(1./Aires));
end
%
if (ordre==2)
  Mpfile=sprintf('MpNC%i_%s',oP,filename);
  Mp=loadMat(Mpfile,Np,Np);
  invMpfile=sprintf('invMpNC%i_%s',oP,filename);
  invMp=loadMat(invMpfile,Np,Np);
end
% Matrices de couplage vitesse-pression
BxFile=sprintf('BxNC%i%i_%s',oU,oP,filename);
Bx=loadMat(BxFile,Np,Nu);
%
ByFile=sprintf('ByNC%i%i_%s',oU,oP,filename);
By=loadMat(ByFile,Np,Nu);