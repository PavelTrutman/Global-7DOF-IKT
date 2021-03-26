% Pavel Trutman
% pavel.trutman@cvut.cz

function M = DHT_GP(M, c, s, i)
  [M1, M2] = DH_GP(M, c, s, i);
  M = M1*M2;
end