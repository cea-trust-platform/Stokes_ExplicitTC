%{
/****************************************************************************
* Copyright (c) 2023, CEA
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
% RHSStokesVortexPolTri.m:
%
% Integrales de $F(nu,alpha,beta)\,P(X,Y)$ dans les triangles qui touchent le centre, 
%        $F(\nu,\alpha,\beta)=-\nu(\alpha^2-1)\rho^{\alpha-2}\evec_\theta
%                               +\beta\rho^{\beta-1}\evec_\rho
%        $\rho$ est la distance a (0,0) et $P(X,Y)\in\{1,X,Y,XY,X^2,Y^2\}$.
%
% SYNOPSIS [RHS_UxC,RHS_UyC]=RHSStokesVortexPolTri(nu,alpha,beta,basis,Ndof)
%
% calcul du second membre pour les aretes qui touchent le centre
%
% INPUT  : - nu
%          - alpha : exposant vitesse
%          - beta  : exposant pression
% OUTPUT : - RHS_UxC = (\int_T F(nu,alpha,beta)\cdot\psi_x) dT
%          - RHS_UyC = (\int_T F(nu,alpha,beta)\cdot\psi_y) dT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHS_UxC,RHS_UyC]=RHSStokesVortexPolTri(nu,alpha,beta,basis,NdofU)
%
global Nbpt CoorNeu xy0som
global Aires NumTri TriEdg xy0tri hmi
global Nbedg CoorMil NumEdg EdgNorm EdgTri
global EdgC
global ordre
%
RHS_UxC=zeros(NdofU,1); RHS_UyC=zeros(NdofU,1);
%
IJK=[1,2,3;2,3,1;3,1,2];
%
% Triangles qui touchent le centre
NbtriC=size(xy0tri,1);
% Aretes qui touchent le centre
NbedgC=size(EdgC,1);
% Angles des aretes qui touchent le centre
ThetaC=zeros(Nbedg,1);
for edg=1:NbedgC
  e=EdgC(edg,1);
  for a=1:2
    iloc=NumEdg(e,a);
    if (iloc~=xy0som)
      X0=CoorNeu(iloc,1)-0.5;
      Y0=CoorNeu(iloc,2)-0.5;
      ThetaC(e,1)=atan2(Y0,X0);
    end
  end
end
nth=1000;
cA=-nu*(alpha^2-1)*0.5; cB=beta*0.5;
invA =cA/alpha; invA1=cA/(1+alpha); invA2=cA/(2+alpha);
invB1=cB/(1+beta); invB2=cB/(2+beta);  invB3=cB/(3+beta);
eps=10.^(-10);
for tri=1:NbtriC
  % triangle du coin
  t=xy0tri(tri,1);
  IGLO=NumTri(t,:);CoorNeuT=CoorNeu(IGLO,:);
  CoorNeuT0=CoorNeu(IGLO,:)-0.5;
  % son aire
  aire=Aires(t);
  % Normales aux faces 
  AGLO=TriEdg(t,:);
  EdgNormT=EdgNorm(AGLO,:);
  for iloc=1:2
   if (EdgTri(AGLO(iloc),1)~=t)
    EdgNormT(iloc,:)=-EdgNormT(iloc,:);
   end
  end
  EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
  % le numero local du centre et ceux des sommets opposes
  i0=xy0tri(tri,2); i1=IJK(i0,2); i2=IJK(i0,3);
  % aretes touchant le centre
  a1=TriEdg(t,i2); th1=ThetaC(a1,1); 
  a2=TriEdg(t,i1); th2=ThetaC(a2,1);
  if ((th1>0)&&(th2<0))
    th2+=2*pi;
  end
  theta=linspace(th1,th2,nth+1);
  dth=theta(1,2:nth+1,1)-theta(1,1:nth);
  costh=cos(theta)'; sinth=sin(theta)';
  cos2th=costh.^2;  sin2th=sinth.^2;
  sincosth=costh.*sinth;
  % sommets opposes au centre  
  X1=CoorNeuT0(i1,1); Y1=CoorNeuT0(i1,2);
  X2=CoorNeuT0(i2,1); Y2=CoorNeuT0(i2,2);
  diffX=abs(X1-X2);
  if (diffX>0.1*eps)
    % deno est toujours positif
    deno=(X1-X2)*sinth-(Y1-Y2)*costh;
    rho=(X1*Y2-X2*Y1)./deno;
    xth=rho.*costh;  yth=rho.*sinth;
  else
    rho=X1./costh; xth=X1*ones(nth+1,1); yth=rho.*sinth;
    fprintf('  diffX=%7.2e\n',diffX);
  end
  rhoA=rho.^alpha; rhoA1=rho.^(alpha+1);
  %
  Fxth  = sinth.*rhoA; intFx =-invA*(dth*(Fxth(1:nth,1)+Fxth(2:nth+1,1)));
  Fyth  = costh.*rhoA; intFy = invA*(dth*(Fyth(1:nth,1)+Fyth(2:nth+1,1)));
  FxXth = sincosth.*rhoA1; intFxX =-invA1*(dth*(FxXth(1:nth,1)+FxXth(2:nth+1,1)));
  FyXth = cos2th.*rhoA1;   intFyX = invA1*(dth*(FyXth(1:nth,1)+FyXth(2:nth+1,1)));
  FxYth = sin2th.*rhoA1;   intFxY =-invA1*(dth*(FxYth(1:nth,1)+FxYth(2:nth+1,1)));
  intFyY=-intFxX;
  %
  rhoB1=rho.^(beta+1); rhoB2=rho.^(beta+2);
  %  
  Gxth  = costh.*rhoB1; intGx=invB1*(dth*(Gxth(1:nth,1)+Gxth(2:nth+1,1)));
  Gyth  = sinth.*rhoB1; intGy=invB1*(dth*(Gyth(1:nth,1)+Gyth(2:nth+1,1)));
  GxXth = cos2th.*rhoB2;  intGxX=invB2*(dth*(GxXth(1:nth,1)+GxXth(2:nth+1,1)));
  GyXth = sincosth.*rhoB2;intGyX=invB2*(dth*(GyXth(1:nth,1)+GyXth(2:nth+1,1)));
  intGxY= intGyX;
  GyYth = sin2th.*rhoB2;  intGyY=invB2*(dth*(GyYth(1:nth,1)+GyYth(2:nth+1,1)));
  intF=[intFx,intFy;intFxX,intFyX;intFxY,intFyY];
  intG=[intGx,intGy;intGxX,intGyX;intGxY,intGyY];
  if (ordre==2)
    cos3th=costh.^3;  sin3th=sinth.^3; sincos2th=sinth.*cos2th; cossin2th=costh.*sin2th;
    rhoA2=rho.^(alpha+2);
    FxXYth = cossin2th.*rhoA2; intFxXY =-invA2*(dth*(FxXYth(1:nth,1)+FxXYth(2:nth+1,1)));
    FyXYth = sincos2th.*rhoA2; intFyXY = invA2*(dth*(FyXYth(1:nth,1)+FyXYth(2:nth+1,1)));
    intFxX2=-intFyXY;
    FyX2th = cos3th.*rhoA2; intFyX2= invA2*(dth*(FyX2th(1:nth,1)+FyX2th(2:nth+1,1)));
    FxY2th = sin3th.*rhoA2; intFxY2=-invA2*(dth*(FxY2th(1:nth,1)+FxY2th(2:nth+1,1)));
    intFyY2=-intFxXY;
    %    
    rhoB3=rho.^(beta+3);
    GxXYth = sincos2th.*rhoB3; intGxXY=invB3*(dth*(GxXYth(1:nth,1)+GxXYth(2:nth+1,1)));
    GyXYth = cossin2th.*rhoB3; intGyXY=invB3*(dth*(GyXYth(1:nth,1)+GyXYth(2:nth+1,1)));
    GxX2th = cos3th.*rhoB3;    intGxX2=invB3*(dth*(GxX2th(1:nth,1)+GxX2th(2:nth+1,1)));
    intGyX2=intGxXY;
    intGxY2=intGyXY;
    GyY2th = sin3th.*rhoB3;    intGyY2=invB3*(dth*(GyY2th(1:nth,1)+GyY2th(2:nth+1,1)));
    %
    intF=[intF;intFxXY,intFyXY;intFxX2,intFyX2;intFxY2,intFyY2];
    intG=[intG;intGxXY,intGyXY;intGxX2,intGyX2;intGxY2,intGyY2];
  end
  %
  intRHS=intF+intG;
  %
  CoorMilT0=CoorMil(AGLO,:)-0.5;
  % Coefficients des fonctions de base de Lagrange
  C_p1=zeros(3,3); coeff=1/(2*aire);
  for i=1:3
    C_p1(i,1)=coeff*(CoorMilT0(i,:)*EdgNormT(i,:)');
  end
  C_p1(:,2:3)=-coeff*EdgNormT;
  %
  if (ordre==2)
    C_p2=zeros(6,6);
    C_p2(1:3,1)=C_p1(:,1).*(2*C_p1(:,1)-1);
    C_p2(1:3,2)=C_p1(:,2).*(4*C_p1(:,1)-1);
    C_p2(1:3,3)=C_p1(:,3).*(4*C_p1(:,1)-1);
    C_p2(1:3,4)=4*C_p1(:,2).*C_p1(:,3);
    C_p2(1:3,5)=2*C_p1(:,2).^2;
    C_p2(1:3,6)=2*C_p1(:,3).^2;
    %
    ik=[2,3,1]; ij=[3,1,2];
    C_p2(4:6,1)=4*C_p1(ik,1).*C_p1(ij,1);
    C_p2(4:6,2)=4*(C_p1(ik,1).*C_p1(ij,2)+C_p1(ik,2).*C_p1(ij,1));
    C_p2(4:6,3)=4*(C_p1(ik,1).*C_p1(ij,3)+C_p1(ik,3).*C_p1(ij,1));
    C_p2(4:6,4)=4*(C_p1(ik,2).*C_p1(ij,3)+C_p1(ik,3).*C_p1(ij,2));
    C_p2(4:6,5)=4*C_p1(ik,2).*C_p1(ij,2);
    C_p2(4:6,6)=4*C_p1(ik,3).*C_p1(ij,3);
  end
  %
  if (strcmp(basis,'LG')==1)
    if (ordre==1)
      RHS_UxC(IGLO,1)+=C_p1*intRHS(:,1);
      RHS_UyC(IGLO,1)+=C_p1*intRHS(:,2);
    end
    %
    if (ordre==2)
      UGLO=[IGLO,AGLO+Nbpt];
      RHS_UxC(UGLO,1)+=C_p2*intRHS(:,1);
      RHS_UyC(UGLO,1)+=C_p2*intRHS(:,2);      
    end
  end
  %
  if (strcmp(basis,'CR')==1)
    C_cr=zeros(3,3); 
    C_cr(:,1)=1-2*C_p1(:,1); 
    C_cr(:,2:3)=-2*C_p1(:,2:3);
    RHS_UxC(AGLO,1)+=C_cr*intRHS(:,1);
    RHS_UyC(AGLO,1)+=C_cr*intRHS(:,2);
  end
  %
  if (strcmp(basis,'CRRT')==1)
    Cx=zeros(3,3); Cx(:,1)=-CoorNeuT0(:,1); Cx(:,2)=ones(3,1);
    Cy=zeros(3,3); Cy(:,1)=-CoorNeuT0(:,2); Cy(:,3)=ones(3,1);
    RHS_tmp=coeff*(Cx*intRHS(:,1)+Cy*intRHS(:,2));    
    RHS_UxC(AGLO,1)+=EdgNormT(:,1).*RHS_tmp;
    RHS_UyC(AGLO,1)+=EdgNormT(:,2).*RHS_tmp;
  end
  %
  if (strcmp(basis,'FS')==1)
    ind=[2,3;3,1;1,2];
    C_fs=zeros(7,6);
    C_fs(1:6,:)=C_p2;
    % Fonction de base triangle
    % coefficients 1,X,Y
    C_fs(7,1)=2-3*sum(C_p1(:,1).^2); 
    C_fs(7,2)=-6*sum(C_p1(:,1).*C_p1(:,2));
    C_fs(7,3)=-6*sum(C_p1(:,1).*C_p1(:,3));
    C_fs(7,4)=-6*sum(C_p1(:,2).*C_p1(:,3));
    C_fs(7,5)=-3*sum(C_p1(:,2).^2); 
    C_fs(7,6)=-3*sum(C_p1(:,3).^2); 
    %
    UGLO=[IGLO,AGLO+Nbpt,Nbpt+Nbedg+t];
    RHS_UxC(UGLO,1)+=C_fs*intRHS(:,1);
    RHS_UyC(UGLO,1)+=C_fs*intRHS(:,2);
  end
  %
  if (strcmp(basis,'FSRT')==1)
    [abdX,abdY,matB]=PiRT1_abdXY(aire,CoorNeuT0,EdgNormT);
    Cxx=zeros(7,6); Cxy=zeros(7,6);
    Cyx=zeros(7,6); Cyy=zeros(7,6);
    BL=matB*C_p1(2:3,:);
    for fb=1:7
      Mx=matB*[abdX(1:2,fb)';abdX(3:4,fb)']; bx=abdX(5:6,fb); dx=matB*abdX(7:8,fb);
      Mx23=Mx*C_p1(2:3,:); bxL=bx'*C_p1(2:3,:);
      Cxx(fb,1:3)=Mx23(1,:); Cxx(fb,1)+=dx(1)+bxL(1)*BL(1,1);
      Cxy(fb,1:3)=Mx23(2,:); Cxy(fb,1)+=dx(2)+bxL(1)*BL(2,1);
      Cxx(fb,2)+=bxL(2)*BL(1,1)+bxL(1)*BL(1,2);
      Cxy(fb,2)+=bxL(2)*BL(2,1)+bxL(1)*BL(2,2);
      Cxx(fb,3)+=bxL(3)*BL(1,1)+bxL(1)*BL(1,3);
      Cxy(fb,3)+=bxL(3)*BL(2,1)+bxL(1)*BL(2,3);
      Cxx(fb,4)+=bxL(2)*BL(1,3)+bxL(3)*BL(1,2);
      Cxy(fb,4)+=bxL(2)*BL(2,3)+bxL(3)*BL(2,2);
      Cxx(fb,5)+=bxL(2)*BL(1,2);
      Cxy(fb,5)+=bxL(2)*BL(2,2);
      Cxx(fb,6)+=bxL(3)*BL(1,3);
      Cxy(fb,6)+=bxL(3)*BL(2,3);
      %     
      My=matB*[abdY(1:2,fb)';abdY(3:4,fb)']; by=abdY(5:6,fb); dy=matB*abdY(7:8,fb);
      My23=Mx*C_p1(2:3,:); byL=by'*C_p1(2:3,:);
      Cyx(fb,1:3)=My23(1,:); Cyx(fb,1)+=dy(1)+byL(1)*BL(1,1);
      Cyy(fb,1:3)=My23(2,:); Cyy(fb,1)+=dy(2)+byL(1)*BL(2,1);
      Cyx(fb,2)+=byL(2)*BL(1,1)+byL(1)*BL(1,2);
      Cyy(fb,2)+=byL(2)*BL(2,1)+byL(1)*BL(2,2);
      Cyx(fb,3)+=byL(3)*BL(1,1)+byL(1)*BL(1,3);
      Cyy(fb,3)+=byL(3)*BL(2,1)+byL(1)*BL(2,3);
      Cyx(fb,4)+=byL(2)*BL(1,3)+byL(3)*BL(1,2);
      Cyy(fb,4)+=byL(2)*BL(2,3)+byL(3)*BL(2,2);
      Cyx(fb,5)+=byL(2)*BL(1,2);
      Cyy(fb,5)+=byL(2)*BL(2,2);
      Cyx(fb,6)+=byL(3)*BL(1,3);
      Cyy(fb,6)+=byL(3)*BL(2,3);
    end
    UGLO=[IGLO,AGLO+Nbpt,Nbpt+Nbedg+t];
    RHS_UxC(UGLO)+=Cxx*intRHS(:,1)+Cxy*intRHS(:,2);
    RHS_UyC(UGLO)+=Cyx*intRHS(:,1)+Cyy*intRHS(:,2);
    %
    %   
    % Integration points
    [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT0);
    rho=sqrt(xyp(:,1).^2+xyp(:,2).^2);
    theta=atan2(xyp(:,2),xyp(:,1));
    cth=cos(theta); sth=sin(theta);
    rha=rho.^(alpha-2); rhb=rho.^(beta-1);
    awp=3*aire*wp';
    FxLOC= nu*(alpha^2-1)*rha.*sth+beta*rhb.*cth;
    FyLOC=-nu*(alpha^2-1)*rha.*cth+beta*rhb.*sth;
    %
    [PiRTxx,PiRTxy,PiRTyx,PiRTyy]=PiRT1(aire,CoorNeuT0,EdgNormT,lambda);
    RHS_UxC(UGLO)-(FxLOC'*PiRTxx+FyLOC'*PiRTxy)';
    RHS_UyC(UGLO)-(FxLOC'*PiRTyx+FyLOC'*PiRTyy)';
  end
 end