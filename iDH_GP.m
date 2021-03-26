% Pavel Trutman
% pavel.trutman@cvut.cz

function [M1, M2] = iDH_GP(M, c, s, i)
  M1 = [ c(i) s(i) 0       0;
        -s(i) c(i) 0       0;
            0    0 1 -M.d(i);
            0    0 0       1];
    
  M2 = [1                0               0 -M.a(i);
        0  cos(M.alpha(i)) sin(M.alpha(i))       0;
        0 -sin(M.alpha(i)) cos(M.alpha(i))       0;
        0                0               0       1];
end