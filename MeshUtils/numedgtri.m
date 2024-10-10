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
% numedgtri.m:
% routine qui indique pour une arete les numeros locaux de ses extremites
% dans les triangles auquels elle appartient. 
%  
% SYNOPSIS LocNumEdgTri=numedgtri();
%          
% OUTPUT Numero locaux des extremites de l'arete et des sommets opposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocNumEdgTri=numedgtri()
  
global NumTri NumEdg EdgTri SomOpp 
%
NbEdg=size(EdgTri,1);
LocNumEdgTri=zeros(NbEdg,2,3);
%
for f=1:NbEdg
  %
  ij=NumEdg(f,:);
  T12=EdgTri(f,:);
  for t=1:2
    Tt=T12(t);
    if (Tt>0)
      NumTriT=NumTri(Tt,:);
      sop=SomOpp(f,t);
      for i=1:3
        iloc=NumTriT(i);
        for e=1:2
          if (iloc==ij(e))
            LocNumEdgTri(f,t,e)=i;
            break;
          end
        end
        %
        if (iloc==sop)
          % Milieu pour l'ordre 2
          LocNumEdgTri(f,t,3)=3+i;
        end
      end % for i=1:3
    end
  end % for t=1:2
end % for f=1:Nbedg
  
  ##for e=1:NbEdg
##  T12=EdgTri(e,:);
##  S12=NumEdg(e,:);
##  for t=1:2
##    Tt=T12(t);
##    if (Tt>0)
##      % sommets du triangle
##      St=NumTri(Tt,:);
##      Sopp=SomSopp(e,t);
##      for s=1:2
##        
##        LocNumEdgTri(e,t,s)=lookup(St,s);
##      end
##      LocNumEdgTri(e,t,3)=lookup(St,Sopp);
##    end
##  end
##end