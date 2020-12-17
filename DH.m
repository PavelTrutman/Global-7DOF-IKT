% Pavel Trutman
% pavel.trutman@cvut.cz

function [M1, M2] = DH(M, th, i)
  M1 = [cos(th(i)) -sin(th(i)) 0      0;
        sin(th(i))  cos(th(i)) 0      0;
                 0           0 1 M.d(i);
                 0           0 0      1];
    
  M2 = [1               0                0 M.a(i);
        0 cos(M.alpha(i)) -sin(M.alpha(i))      0;
        0 sin(M.alpha(i))  cos(M.alpha(i))      0;
        0               0                0      1];
end