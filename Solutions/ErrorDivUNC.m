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
% Calcul de l'erreur sur l'approximation de la divergence pour un vecteur P1 ou P2 NC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Eud=ErrorDivUNC(Uxh,Uyh)
  global Nbedg EdgTri LgEdg2 EdgNorm NumEdg SomOpp RefEdg
  global TriEdg NumTri
  %
  Eud=0;
  %
  IJK=[1,2,3;2,3,1;3,2,1];
  EdgI=find(RefEdg==0);
  NbedgI=size(EdgI,1);
  % on limite le calcul aux aretes interieures
  for fi=1:NbedgI
    f=EdgI(fi);
    ILOCF=NumEdg(f,:); hF2=LgEdg2(f,1);
    s1=ILOCF(1); s2=ILOCF(2);
    %    
    EdgNormF=EdgNorm(f,:)';
    %
    TT=EdgTri(f,:);
    %
    t1=TT(1); ALOC1=TriEdg(t1,:); Sf1=SomOpp(f,1);
    ILOC1=NumTri(t1,:);
    i1=find(ILOC1==Sf1); j1=IJK(i1,2); k1=IJK(i1,3);
    js1=find(ILOCF==ILOC1(j1)); ks1=find(ILOCF==ILOC1(k1));
    Uxh1=Uxh(ALOC1); Uyh1=Uyh(ALOC1);
    UhF1=zeros(3,2);
    UhF1(js1,1)=-Uxh1(j1)+Uxh1(i1)+Uxh1(k1); UhF1(js1,2)=-Uyh1(j1)+Uyh1(i1)+Uyh1(k1);
    UhF1(ks1,1)=-Uxh1(k1)+Uxh1(i1)+Uxh1(j1); UhF1(ks1,2)=-Uyh1(k1)+Uyh1(i1)+Uyh1(j1);
    UhF1(3,1)=Uxh1(i1); UhF1(3,2)=Uyh1(i1);
    %
    t2=TT(2); ALOC2=TriEdg(t2,:); Sf2=SomOpp(f,2);
    ILOC2=NumTri(t2,:);
    i2=find(ILOC2==Sf2); j2=IJK(i2,2); k2=IJK(i2,3);
    js2=find(ILOCF==ILOC1(j2)); ks2=find(ILOCF==ILOC2(k2));
    Uxh2=Uxh(ALOC2); Uyh2=Uyh(ALOC2);
    UhF2=zeros(3,2);
    UhF2(js2,1)=-Uxh2(j2)+Uxh2(i2)+Uxh2(k2); UhF2(js2,2)=-Uyh2(j2)+Uyh2(i2)+Uyh2(k2);
    UhF2(ks2,1)=-Uxh2(k2)+Uxh2(i2)+Uxh2(j2); UhF2(ks2,2)=-Uyh2(k2)+Uyh2(i2)+Uyh2(j2);
    UhF2(3,1)=Uxh2(i2); UhF2=Uyh2(i2);
    %
    wp=zeros(1,3); wp(1,1)=1/6; wp(1,2)=wp(1); wp(1,3)=1-wp(1)-wp(2);
    Eud=Eud+wp*(((UhF1-UhF2)*EdgNormF).^2)/hF2;
      %
  end
  Eud=sqrt(Eud);