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
%  nu*A*U(n+1)+Bt*P(n+1) = F-nu*Ab*Ub
% -la*B*U(n+1)+(1/nu)*M*P(n+1)=(1/nu)*M*P(n)+la*Bb*Ub
%
% ALGO:
% (1) nu*(A+la*Bt*invM*B)*U(n+1)=F-nu*(Ab+la*Bt*invM*Bb)*Ub-Bt*P(n) = RHS_U
% (2) P(n+1)=P(n)+nu*la*invM*(Bb*Ub+B*U(n+1))
%
% SolveStokesExplicitTC_it
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxh,Uyh,Ph,Eud,nit]=SolveStokesExplicitTC_it(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,Ux0,Uy0,DofU,DofUB,Ph0,dPh0,Vol)
  global nu lambda nmax pas epsTC
  %
  invLa=1/lambda;
  invNu=1/nu; invNu2=invNu^2;
  %
  Nuxy=size(RHS_Ux,1);
  Uxh=zeros(Nuxy,1); Uyh=Uxh; AUx=Uxh; AUy=Uyh;
  UxhB=Ux0(DofUB,1); UyhB=Uy0(DofUB,1);
  Uxh(DofUB,1)=UxhB;
  Uyh(DofUB,1)=UyhB;
  %
  Nu=size(DofU,1); Np=size(Ph0,1);
  %
  Kubb=Ku(DofUB,DofUB);
  normU1B2=UxhB'*Kubb*UxhB+UyhB'*Kubb*UyhB;
  Ku0=Ku(DofU,DofU); Zu=sparse(Nu,Nu); Bx0=Bx(:,DofU); By0=By(:,DofU); 
  %
  Kub=Ku(DofU,DofUB); Bxb=Bx(:,DofUB); Byb=By(:,DofUB); 
  nuUxB=nu*UxhB; nuUyB=nu*UyhB;
  %
  nuAUxb=Kub*nuUxB; nuAUyb=Kub*nuUyB;
  %
  Pb=lambda*invMp*(Bxb*nuUxB+Byb*nuUyB);
  %
  RHS_Ux0=RHS_Ux(DofU,1)-nuAUxb;
  RHS_Uy0=RHS_Uy(DofU,1)-nuAUyb;
  % 
  Kmat=[Ku0,Zu;Zu,Ku0]; Bmat=[Bx0,By0];
  invMpaB=lambda*invMp*Bmat;
  %
  Amat=Kmat+Bmat'*invMpaB; 
  Ph=Ph0;
  %
  normDivU=invLa*invNu*sqrt(dPh0'*Mp*dPh0);
  %
  nit=1; Sol=[nu*Ux0(DofU,1);nu*Uy0(DofU,1)];
  %
  stop=0;
  Bx0Ph=Bx0'*Ph; By0Ph=By0'*Ph;
  %  
  RHS_UxIt=RHS_Ux0-Bx0'*Pb;
  RHS_UyIt=RHS_Uy0-By0'*Pb;
  %
  while ( (nit<nmax) && (stop==0) )
    nit+=1; normDivU0=normDivU; Ph0=Ph; Sol0=Sol;
    RHS_UxSol=RHS_UxIt-Bx0Ph;
    RHS_UySol=RHS_UyIt-By0Ph;
    RHS=[RHS_UxSol;RHS_UySol];
    Sol=Amat\RHS; % Attention : Sol=\nu[Uxh;Uyh]_\omega;
    Ph=Ph0+invMpaB*Sol+Pb;
    %
    moyPh=UnP'*Mp*Ph/Vol;
    Ph=Ph-moyPh;
    dP=Ph-Ph0;
    normDivU=invLa*invNu*sqrt(dP'*Mp*dP);
    ratio=normDivU/normDivU0;
    if (ratio>1)    
      fprintf('Iteration %i: ||div u_h^n||_0/||div u_h^{n-1}||_0=%7.2e\n',nit,ratio); 
    end
    %
    Uxit=Sol(1:Nu,1); Uyit=Sol(Nu+1:2*Nu,1);
    %   
    Bx0Ph=Bx0'*Ph; RHS_Ux0It=RHS_Ux0-Bx0Ph;
    By0Ph=By0'*Ph; RHS_Uy0It=RHS_Uy0-By0Ph;
    %
    normH1U=normU1B2+invNu2*(RHS_Ux0It'*Uxit+RHS_Uy0It'*Uyit)+2*invNu2*(nuAUxb'*Uxit+nuAUyb'*Uyit);
    normH1U=sqrt(normH1U);
    % VERIFICATION !!
    %Uxh(DofU,1)=invNu*Uxit; Uyh(DofU,1)=invNu*Uyit;
    %[normDivUh,normRotUh,normH1Uh]=normDivRotU(Uxh,Uyh);
    %
    if (normDivU<(epsTC*normH1U))
      fprintf('Iteration %i: ||div u_h^n||_0/|u_h^n|_1=%7.2e\n',nit,normDivU/normH1U);
      stop=1;
    end
    %
  end
  %
  Eud=normDivU;
  Uxh(DofU,1)=invNu*Sol(1:Nu,1); 
  Uyh(DofU,1)=invNu*Sol(Nu+1:2*Nu,1);
  fprintf('ExplicitTC: lambda=%7.2e, ||div Uh||_0/|Uh|_1=%7.2e\n',lambda,Eud/normH1U);