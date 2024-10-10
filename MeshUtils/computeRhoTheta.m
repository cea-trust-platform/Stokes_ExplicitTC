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
% computeRhoTheta.m:
%
% Calcul des coordonnées polaires
%
% SYNOPSIS [rho2,rho,theta,xy0tri,xy0edg,xy0som]=computeRhoTheta(xy0,ordre)
%
% INPUT  : - xy0    : coordonnees du centre du vortex
% OUTPUT : - rho2   : distance des points a xy0 au carre
%          - rho    : distance des points a xy0
%          - theta  : angle polaire dans le repere (xy0,X,Y)
%          - xy0tri : triangle(s) dans lequel se trouve xy0
%          - xy0edg : arete sur laquelle se trouve xy0
%          - xy0som : sommet sur lequel se trouve xy0
%          - som    : sommet le plus proche de xy0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho2,rho,theta,xy0tri,xy0edg,xy0som,som]=computeRhoTheta(xy0)
  global meshstep ordre
  global Nbtri NumTri TriEdg Aires
  global Nbpt CoorNeu
  global Nbedg NumEdg EdgNorm EdgTri CoorMil
  %
  xy0tri=[]; xy0edg=[]; xy0som=0;
  %
  XY=CoorNeu;
  if (ordre==2)
    XY=[XY;CoorMil];
  end
  %
  Ndof=size(XY,1);
  Gom=xy0.*ones(Ndof,2);
  XY0=XY-Gom;
  % Calcul de rho
  rho2=XY0(:,1).^2+XY0(:,2).^2;
  rho=sqrt(rho2);
  % Localisation de xy0 dans le maillage
  % Recherche pour voir si un noeud a ces coordonnees
  eps=1.e-15;
  % Indice du sommet le plus proche de xy0
  som=0;
  diff=1;
  for i=1:Ndof
    diffi=rho(i)/meshstep;
    if (diffi<=eps)
      if (i<=Nbpt)
        xy0som=i; 
        som=xy0som;
      else
        % le degre de liberte i
        %  est un numero d'arete
        xy0edg=i-Nbpt;
      end
      break;
    end
    if (diffi<=diff)
      som=i;
      diff=diffi;
    end
  end
  %fprintf('Distance min = %7.2e, diff = %7.2e\n',rho(som),diff);
  %
  % Calcul de theta
  theta=atan2(XY0(:,2),XY0(:,1));
  %
  % Pour chaque sommet, numeros de ses triangles
  TriSom=[];
  NtriSom=zeros(Nbpt,1);
  for t=1:Nbtri
    for i=1:3
      vi=NumTri(t,i);
      Nti=NtriSom(vi)+1;
      TriSom(vi,Nti)=t;
      NtriSom(vi)=Nti;
    end
  end
  %
  % Pour chaque sommet, numeros de ses aretes
  EdgSom=[];
  NedgSom=zeros(Nbpt,1);
  for e=1:Nbedg
    for i=1:2
      vi=NumEdg(e,i);
      Nei=NedgSom(vi)+1;
      EdgSom(vi,Nei)=e;
      NedgSom(vi)=Nei;
    end
  end
  %
  if (xy0edg~=0)
    % Indices des triangles ayant xy0edg pour face
    EdgTrixy0=EdgTri(xy0edg,:)';
    Nt=size(EdgTrixy0,1);
    xy0tri=zeros(Nt,2);
    xy0tri(:,1)=EdgTrixy0;
  end
  %
  if (xy0som~=0)
    % Indices des triangles ayant xy0 pour sommet
    Nt=NtriSom(xy0som,1);
    xy0tri=zeros(Nt,2);
    xy0tri(:,1)=TriSom(xy0som,1:Nt)';
    for t=1:Nt
      TriT=xy0tri(t,1);
      % indice local de xy0 dans le triangle
      IGLO=NumTri(TriT,:);
      iloc=0;
      for i=1:3
        if (IGLO(i)==xy0som)
          iloc=i;
          xy0tri(t,2)=iloc;
          break
        end
      end
    end
    % Indices des aretes ayant xy0 pour extremite
    Ne=NedgSom(xy0som,1);
    xy0edg=zeros(Ne,2);
    xy0edg(:,1)=EdgSom(xy0som,1:Ne)';
    for e=1:Ne
      EdgT=xy0edg(e,1);
      % indice local de xy0 dans l'arete
      AGLO=NumEdg(EdgT,:);
      aloc=0;
      for i=1:2
        if (AGLO(i)==xy0som)
          aloc=i;
          xy0edg(e,2)=aloc;
          break
        end
      end
    end
  end
  %
  % Indices des sommets proches de xy0
  if (xy0som==0)
    % Nombre de triangles contenant som
    Nt=NtriSom(som,1);
    Tri=TriSom(som,1:Nt);
    %
    for t=1:Nt
      TriT=Tri(1,t);
      IGLO=NumTri(TriT,:);
      CoorNeuT=CoorNeu(IGLO,:)
    end
    for t=1:Nt
      TriT=Tri(1,t);
      IGLO=NumTri(TriT,:);
      NumtriS(t,:)=IGLO;
      AGLO=TriEdg(TriT,:);
      CoorNeuT=CoorNeu(IGLO,:);
      EdgNormT=EdgNorm(AGLO,:);
      aire=Aires(TriT);
      for iloc=1:2
        if (EdgTri(AGLO(iloc),1)~=TriT)
          EdgNormT(iloc,:)=-EdgNormT(iloc,:);
        end
      end
      EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
      CoorMilT=CoorMil(AGLO,:);
      % Calcul des coordonnees barycentriques
      coeff=0.5/Aires(TriT);
      Mat=coeff*[-EdgNormT,sum(CoorMilT.*EdgNormT,2)];
      lambda=[xy0,1]*Mat';
      lneg=0;
      for i=1:3
        if (lambda(i)<0)
          % xy0 n'est pas dans TriT
          lneg=1;
          break;
        end
      end
      % on a trouve le triangle
      if (lneg==0)
        lb=0;
        iloc=0;
        for i=1:3
          if (lambda(i)>lb)
            iloc=i;
            lb=lambda(i);
          end
        end
        xy0tri=[TriT,iloc];
        % le point xy0 est-il sur une arete ?
        for i=1:3
          if (lambda(i)<eps)
            edg=ALGO(i);
            aloc=1;
            if (NumEdg(edg,2)==som)
              aloc=2;
            end
            xy0edg=[edg,aloc];
            break;
          end
        end              
        break;
      end % end_if (lneg==0)
    end % end_for t=1:Nt
  end % end_if (xy0som==0)