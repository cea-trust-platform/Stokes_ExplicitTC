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
% loadMatStokesSV.m:
%
% Lecture des matrices de masse, de raideur et de couplage
%             EF P1-P0 (ordre=1) ou P2-P1dg (ordre=2) 
%
% SYNOPSIS [Ku,Mu,Mp,invMp,Bx,By] = loadMatStokesSV(filename)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ku,Mu,Mp,invMp,Bx,By] = loadMatStokesSV(filename)
  global ordre Nbpt Nbtri Nbedg
  Nu=Nbpt; Np=Nbtri;
  oU=ordre; oP=ordre-1;
  if (ordre==2)
    Nu+=Nbedg;Np*=3;
  end
  Kufile=sprintf('KuSV%i_%s',oU,filename);
  Ku=loadMat(Kufile,Nu,Nu);
  Mufile=sprintf('MuSV%i_%s',oU,filename);
  Mu=loadMat(Mufile,Nu,Nu);
  Mpfile=sprintf('MpSV%i_%s',oP,filename);
  Mp=loadMat(Mpfile,Np,Np);
  invMpfile=sprintf('invMpSV%i_%s',oP,filename);
  invMp=loadMat(invMpfile,Np,Np);
  Bxfile=sprintf('BxSV%i%i_%s',oU,oP,filename);
  Bx=loadMat(Bxfile,Np,Nu);
  Byfile=sprintf('BySV%i%i_%s',oU,oP,filename);
  By=loadMat(Byfile,Np,Nu);