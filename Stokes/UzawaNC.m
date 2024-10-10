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
% UzawaNC.m:
%
% Algorithme de resolution du systeme mixte 
% ( nu*A     0   Bx^T)(Uxh) =(RHS_Ux)
% (    0  nu*A   By^T)(Uyh) =(RHS_Uy)
% (  -Bx   -By     0 )(Ph ) =(0     )
%
% SYNOPSIS [Uxh,Uyh,Ph]=UzawaNC(Ku0,Bx0,By0,Mp,invMp,UnP,RHS_Ux0,RHS_Uy0,UxB,UyB,DofU,DofUB,Vol)
%          
% Uzawa : on resout B(A^{-1})B^T\,Ph=nu*RHS_p+B(A^{-1})RHS_U par un GCP, Cholesky pour inverser A
%
% INPUT : - Ku0(ndofU0,ndofU0): matrice de raideur interne (une composante)
%         - Bx0(ndofP,ndofU0) : matrice de couplage vitesse x-pression
%         - By0(ndofP,ndofU0) : matrice de couplage vitesse y-pression
%         - Mp(ndofP,ndofP) : matrice de masse pression
%         - invMp(ndofP,1)  : inverse matrice de masse pression
%         - UnP(ndofP,1)    : representation discrete de la fonction 1 sur \Omega.
%         - RHS_Ux0(ndofU0,1) : second membre vitesse x
%         - RHS_Uy0(ndofU0,1) : second membre vitesse y
%         - DofU(:,1)    : numeros des degres de liberte non elimines pour Ux et Uy
%         - Vol          : aire totale
% OUTPUT: - Uxh0(ndofU0,1) : solution vitesse x
%         - Uyh0(ndofU0,1) : solution vitesse y
%         - Ph(ndofP,1)  : solution pression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxh0,Uyh0,Ph,nit]=UzawaNC(Ku0,Bx0,By0,Mp,invMp,UnP,RHS_Ux0,RHS_Uy0,UxB,UyB,DofU,DofUB,Vol)
  global nu eps nitMAX mi ordre
  %
  ndofU0=size(RHS_Ux0,1); ndofU=size(DofU,1);
  Uxh0=zeros(ndofU0,1); Uyh0=Uxh0;
  Uxh=zeros(ndofU,1); Uyh=Uxh;
  %
  Ku=Ku0(DofU,DofU);
  RHS_Ux=RHS_Ux0(DofU,1)-nu*Ku0(DofU,DofUB)*UxB(DofUB);
  RHS_Uy=RHS_Uy0(DofU,1)-nu*Ku0(DofU,DofUB)*UyB(DofUB); 
  RHS_p=nu*(Bx0(:,DofUB)*UxB(DofUB)+By0(:,DofUB)*UyB(DofUB));
  Bx=Bx0(:,DofU); By=By0(:,DofU);
  %%%%%%%%%%%%%%%%
  % Renumerotation 
  %%%%%%%%%%%%%%%%
  sU = symrcm(Ku);
  Ku_s=Ku(sU,sU);
  %%%%%%%%%%%
  % Cholesky 
  %%%%%%%%%%%
  Uchol=chol(Ku_s);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Renumerotation "espace-composante"
  %             en "composante-espace"
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ndofVit=2*ndofU;
  Projx=sparse(ndofU,ndofVit);
  Projy=sparse(ndofU,ndofVit);
  for i=1:ndofU
    Projx(i,2*i-1)=1;
    Projy(i,2*i  )=1;
  end 
  Rmat=Projx'*Uchol*Projx+Projy'*Uchol*Projy;
  Bmat=Bx(:,sU)*Projx+By(:,sU)*Projy;
  Fvec=Projx'*RHS_Ux(sU)+Projy'*RHS_Uy(sU);
  %
  % Vitesse initiale
  Uvec = Rmat\(Rmat'\Fvec);
  Div=Bmat*Uvec; 
  Err0=sqrt(Div'*invMp*Div);
  % Second membre pour le systeme de resolution de la pression
  %
  FPvec = RHS_p+Bmat*Uvec;
  % Pression initiale
  Ph=invMp*FPvec;
  % Projection sur pressions a vitesses nulles
  Ph=Ph-(UnP'*FPvec/Vol);
  % Calcul du residu initial
  Ap=BKBtp(Rmat,Bmat,Ph);
  r=FPvec-Ap;
  %
  % direction de descente initiale
  z=invMp*r;
  d=z;
  d-=(UnP'*r/Vol); % d=d-(UnP'*Mp*z/Vol);
  % norme du residu au carre
  rr=r'*r;
  %
  rz=r'*z;
  %
  %eps=1.e-14;
  eps2=eps^2*rr;
  nit=0;
  while ((rr>eps2)&&(nit<nitMAX))
    Ad=BKBtp(Rmat,Bmat,d);
    alPha=rz/(Ad'*d);
    Ph+=alPha*d;
    r-=alPha*Ad;
    z=invMp*r;
    rz1=r'*z;
    beta=rz1/rz;
    d=z+beta*d;
    d-=(UnP'*r/Vol);% d=d-(UnP'*Mp*d/Vol);
    rr=r'*r;
    rz=rz1;
    nit+=1;
  end
  Ph=Ph-(UnP'*Mp*Ph/Vol);
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Resolution de la vitesse
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  Uvec=Uvec-Rmat\(Rmat'\(Bmat'*Ph));
  Uvec=(1/nu)*Uvec;
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Renumerotations inverses
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  indU=linspace(1,ndofU,ndofU)';
  Cx=2*indU-1; Cy=2*indU;
  Uxh(sU)=Uvec(Cx); 
  Uyh(sU)=Uvec(Cy);
  Uxh0(DofU)=Uxh; Uxh0(DofUB)=UxB(DofUB);
  Uyh0(DofU)=Uyh; Uyh0(DofUB)=UyB(DofUB);
  %
end
