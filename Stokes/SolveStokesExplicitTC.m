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
%  nu*A*U+Bt*P = F-nu*Ab*Ub
% -la*B*U+(1/nu)*M*P=(1/nu)*M*P0+la*Bb*Ub
%
% ALGO:
% (1) nu*(A+la*Bt*invM*B)*U=F-nu*(Ab+la*Bt*invM*Bb)*Ub-Bt*P0 = RHS_U
% (2) P=P0+nu*la*invM*(Bb*Ub+B*U)
%
% SolveStokesExplicitTC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxh,Uyh,Ph,Eud]=SolveStokesExplicitTC(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,Uxh0,Uyh0,DofU,DofUB,Ph0,Vol)
  global nu lambda
  %
  invLa=1/lambda;
  invNu=1/nu; invNu2=invNu^2;
  %
  Nuxy=size(RHS_Ux,1);
  Uxh=zeros(Nuxy,1); Uyh=Uxh;
  %
  Nu=size(DofU,1); Np=size(Ph0,1);
  %  
  UxhB=Uxh0(DofUB,1); UyhB=Uyh0(DofUB,1);
  Uxh(DofUB,1)=UxhB;
  Uyh(DofUB,1)=UyhB;
  %
  Kubb=Ku(DofUB,DofUB);
  normU1B2=UxhB'*Kubb*UxhB+UyhB'*Kubb*UyhB;
  Ku0=Ku(DofU,DofU); Zu=sparse(Nu,Nu); Bx0=Bx(:,DofU); By0=By(:,DofU); 
  %
  nuUxB=nu*UxhB; nuUyB=nu*UyhB;
  Kub=Ku(DofU,DofUB); Bxb=Bx(:,DofUB); Byb=By(:,DofUB);
  Pb=lambda*invMp*(Bxb*nuUxB+Byb*nuUyB);
  %
  nuAUxb=Kub*nuUxB; nuAUyb=Kub*nuUyB;
  %
  RHS_Ux0=RHS_Ux(DofU,1)-nuAUxb-Bx0'*Pb; 
  RHS_Uy0=RHS_Uy(DofU,1)-nuAUyb-By0'*Pb; 
  %
  Kmat=[Ku0,Zu;Zu,Ku0]; Bmat=[Bx0,By0];
  invMpaB=lambda*invMp*Bmat;
  %
  Amat=Kmat+Bmat'*invMpaB; 
  RHS_Ux0=RHS_Ux0-Bx0'*Ph0;
  RHS_Uy0=RHS_Uy0-By0'*Ph0;
  RHS=[RHS_Ux0;RHS_Uy0];
  Sol=Amat\RHS;
  Uxh0=Sol(1:Nu,1); Uyh0=Sol(Nu+1:2*Nu,1);
  Ph=Ph0+invMpaB*Sol+Pb;
  %
  moyPh=UnP'*Mp*Ph/Vol;
  Ph=Ph-moyPh;
  %
  Uxh(DofU,1)=invNu*Uxh0;
  Uyh(DofU,1)=invNu*Uyh0;
  %
  dP=Ph-Ph0;
  %
  normH1U=normU1B2+invNu2*(RHS_Ux0'*Uxh0+RHS_Uy0'*Uyh0)+2*invNu2*(nuAUxb'*Uxh0+nuAUyb'*Uyh0);
  normH1U=sqrt(normH1U);
  %
  Eud=invNu*invLa*sqrt(dP'*Mp*dP);
  fprintf('ExplicitTC: lambda=%7.2e, ||div Uh||_0/|Uh|_1=%7.2e\n',lambda,Eud/normH1U);