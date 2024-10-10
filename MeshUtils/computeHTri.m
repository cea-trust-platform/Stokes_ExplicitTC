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
% SYNOPSIS [DiaTri,invDiaTri,LgEdg,invLgEdg,sigTri]=computeHTri( )
%          Formule des sinus : h_T=a*b*c/2|T|
%  
% OUTPUT - DiaTri(Nbtri,1)  : diametre des triangles
%        - invDiaTri(Nbtri,1)  : inverse des diametres des triangles
%        - LgEdg(Nbedg,1) : longueur des aretes
%        - invLgEdg(Nbedg,1) : inverse des longueurs des aretes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DiaTri,invDiaTri,LgEdg,invLgEdg,sigTri]=computeHTri()
%
global Nbtri TriEdg Aires
global LgEdg2
%
invDiaTri=zeros(Nbtri,1); DiaTri=zeros(Nbtri,1);
PerTri=zeros(Nbtri,1); sigTri=zeros(Nbtri,1);
LgEdg=sqrt(LgEdg2);
invLgEdg=1./LgEdg;
%
for t=1:Nbtri
  DiaTri(t)=1./(2*Aires(t));
  for i=1:3
      a=TriEdg(t,i);
      DiaTri(t)=DiaTri(t)*LgEdg(a);
      PerTri(t)+=LgEdg(a);
  end
end
invDiaTri=1./DiaTri;
rhoTri=4*Aires./PerTri;
sigTri=DiaTri./rhoTri;
%