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
% DivDivxyLGTriangle.m:
%
% Matrices de la divergence EF de Lagrange P1 et P2 pour un triangle donne
%
% P1- : \int\div\vvec\div\uvec
% (i,j)=(\pa_x\lambda_i,\pa_x\lambda_j)_{0,T}
% Dyy(i,j)=(\pa_y\lambda_i,\pa_y\lambda_j)_{0,T}
% Dxy(i,j)=(\pa_x\lambda_i,\pa_y\lambda_j)_{0,T}
% Dyx(i,j)=(\pa_y\lambda_i,\pa_x\lambda_j)_{0,T}
% 
% Pour tester : aire=0.5; EdgNormT=[1,1;-1,0;0,-1];
%
% SYNOPSIS [Dxx,Dyy,Dxy,Dyx]=DivDivxyLGTriangle(aire,EdgNormT)
%          
% INPUT  - aire          : aire du triangle
%        - EdgNormT(3,2) : coordonnees (Nx, Ny) des vecteurs "face-normale" opposes aux sommets locaux
% OUTPUT - Dxx,Dyy,Dxy,Dyx : matrices de couplage vitesse-pression locales pour EF de Lagrange P1 P2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dxx,Dyy,Dxy,Dyx]=DivDivxyLGTriangle(aire,EdgNormT)
  global ordre
  if (ordre==1)
    coeff=0.25/aire;
    EN=coeff*EdgNormT;
    Dxx=EN(:,1)*EdgNormT(:,1)';
    Dxy=EN(:,1)*EdgNormT(:,2)';
    Dyx=EN(:,2)*EdgNormT(:,1)';
    Dyy=EN(:,2)*EdgNormT(:,2)';
  endif
  %
  if (ordre==2)
    
  endif
