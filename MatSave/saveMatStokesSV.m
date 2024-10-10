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
% saveMatStokesSV.m:
%
% Sauvebgarde des matrices du probleme de Stokes, discretisation SV
%
% SYNOPSIS saveMatStokesSV(filename)
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
function saveMatStokesSV(filename)
%
% VARIABLES GLOBALES
%
global Nbpt CoorNeu
global Nbtri Aires NumTri TriEdg
global Nbedg LgEdg2 EdgNorm EdgTri
global ordre
%
for ordre=1:2
  fprintf(' Fabrication des matrices SV ordre=%i.\n',ordre);
   [Ku,Mu,Mp,invMp,Bx,By] = MatStokesSV();
   fprintf(' Sauvegarde des matrices SV ordre=%i.\n',ordre);
   fileMuSV =sprintf('MuSV%i_%s',ordre,filename);
   saveMat(Mu,fileMuSV);
   %
   fileKuSV =sprintf('KuSV%i_%s',ordre,filename);
   saveMat(Ku,fileKuSV);
   %
   fileMpSV =sprintf('MpSV%i_%s',ordre-1,filename);
   saveMat(Mp,fileMpSV);
   %
   fileInvMpSV =sprintf('invMpSV%i_%s',ordre-1,filename);
   saveMat(invMp,fileInvMpSV);
   %
   fileBxSV =sprintf('BxSV%i%i_%s',ordre,ordre-1,filename);
   saveMat(Bx,fileBxSV);
   %
   fileBySV =sprintf('BySV%i%i_%s',ordre,ordre-1,filename);
   saveMat(By,fileBySV);
   %
   Sp=MatStabDGLG(ordre-1);
   fileSpSV =sprintf('SpSV%i%i_%s',ordre,ordre-1,filename);
   saveMat(Sp,fileSpSV);
end