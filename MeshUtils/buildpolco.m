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
% Calcul et stockages des coordonnees polaires centrees en xy0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
global Nbpt CoorNeu CoorNeu2 RefNeu
global Nbtri CoorBary Aires NumTri NumTri2 TriEdg
global Nbedg CoorMil RefEdg LgEdg2 EdgNorm EdgTri NumEdg
global ordre meshstep
%
ordre=2;
meshsteps=[0.1,0.05,0.025,0.0125,0.00625];
global xy0=[0.5,0.5];
global meshstep
%
for mi=1:5
  meshstep=meshsteps(mi);
  filename =sprintf('Vtx_square%i',mi);
  %buildmeshfiles(filename);
  [CoorNeu,CoorNeu2,CoorBary,RefNeu,RefNeu2,NumTri,NumTri2,RefTri,RefTri2,NumEdg,NumEdgB,CoorMil,...
         RefEdg,RefEdgB,TriEdg,EdgTri,SomOpp,LgEdg2,EdgNorm,Aires]=readmeshfiles(filename);
  %
  Nbpt = size(CoorNeu,1); 
  Nbtri= size(NumTri,1) ; 
  Nbedg= size(NumEdg,1) ;
  [rho2,rho,theta,xy0tri,xy0edg,xy0som,som]=computeRhoTheta(xy0); 
  writepolco(rho2,rho,theta,xy0tri,xy0edg,xy0som,filename);
end