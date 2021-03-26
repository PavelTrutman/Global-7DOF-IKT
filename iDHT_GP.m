% Pavel Trutman
% pavel.trutman@cvut.cz

function M = iDHT_GP(M, c, s, i)
  [M1, M2] = iDH_GP(M, c, s, i);
  M = M2*M1;
end