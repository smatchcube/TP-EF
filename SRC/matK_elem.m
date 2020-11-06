function [Kel] = matK_elem(S1, S2, S3, zone)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%         zone : 1 ou 2 en fonction de la zone du triangle concerné 
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) Utilisation d une quadrature a 3 point d ordre 2
%    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;


% calcul de la matrice de raideur
% -------------------------------
B = [x2-x1, x3-x1;
     y2-y1, y3-y1];    
     
if zone == 1
  Kel = ones(3,3) * (1/6) * ( sigma_1(num2cell(B * [1/6; 1/6] + [x1; y1]){:})
                            + sigma_1(num2cell(B * [2/3; 1/6] + [x1; y1]){:})
                            + sigma_1(num2cell(B * [1/6; 2/3] + [x1; y1]){:}));
else 
  Kel = ones(3,3) * (1/6) * ( sigma_2(num2cell(B * [1/6; 1/6] + [x1; y1]){:})
                            + sigma_2(num2cell(B * [2/3; 1/6] + [x1; y1]){:})
                            + sigma_2(num2cell(B * [1/6; 2/3] + [x1; y1]){:}));
end;
                               
diagrent_w = [-1, 1, 0; 
              -1, 0, 1];

for i=1:3
  for j=1:3
    Kel(i,j) *= abs(D) * dot(inverse(transpose(B))*diagrent_w(:,i), 
                             inverse(transpose(B))*diagrent_w(:,j));
  end; % j
end; % i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020
