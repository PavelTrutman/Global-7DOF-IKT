tic:=proc() time(); end proc:
toc:=proc(t) time()-t; end proc:

M2L:=proc(M) 
  convert(convert(M,Vector),list);
end proc:

dhTs := proc(i)
  local M1, M2;
  M1:=Matrix(4,4,[[+cat(`c`,i), -cat(`s`,i), 0,          0],
                 [ +cat(`s`,i), +cat(`c`,i), 0,          0],
                 [           0,           0, 1, cat(`d`,i)],
                 [           0,           0, 0,          1]]);
  M2:=Matrix(4,4,[[1,        0,           0,cat(`a`,i)],  
                 [ 0,+cat(`lambda`,i),-cat(`mu`,i),         0],
                 [ 0,+cat(`mu`,i),+cat(`lambda`,i),         0],
                 [ 0,          0,          0,         1]]);
  [M1,M2];
end proc:

dhInvs := proc(M)
  local M1, M2;
  M1 := M[1];
  M2 := M[2];
  [simplify(MatrixInverse(M2),{M2[3,2]^2+M2[3,3]^2=1}),
   simplify(MatrixInverse(M1),{M1[1,1]^2+M1[2,1]^2=1})];
end proc:

dhT := proc(i)
local M;
  M:=dhTs(i);
  M[1].M[2];
end proc:

dhInv := proc(M)
  simplify(MatrixInverse(M),{M[1,1]^2+M[2,1]^2=1,M[3,2]^2+M[3,3]^2=1});
end proc:

q2R := proc(u)
local R;
  R := <<u(1)^2+u(2)^2-u(3)^2-u(4)^2 | 2*(u(2)*u(3)-u(1)*u(4)) | 2*(u(2)*u(4)+u(1)*u(3))>, <2*(u(2)*u(3)+u(1)*u(4)) | u(1)^2-u(2)^2+u(3)^2-u(4)^2 | 2*(u(3)*u(4)-u(1)*u(2))>, <2*(u(2)*u(4)-u(1)*u(3)) | 2*(u(3)*u(4)+u(1)*u(2)) | u(1)^2-u(2)^2-u(3)^2+u(4)^2>>:
  R := R/(u(1)^2 + u(2)^2 + u(3)^2 + u(4)^2);
end proc:

R2q:=proc(R)
  local q;
  q := <Trace(R)+1, R(3, 2)-R(2, 3), R(1, 3)-R(3, 1), R(2, 1)-R(1, 2)>/(2*sqrt(Trace(R)+1));
  q := q/norm(q, 2);
end proc:

CS2CSR:=proc(th,e)
  local c,s;
  local t, te;
  t := tan(th/2);
  te := round(evalf(t*e))/e;
  c := -(te^2 - 1)/(te^2 + 1);
  s := 2*te/(te^2 + 1);
  c, s;
end proc:

E2ER:=proc(E,e)
  local R, t;
  local te, q, qe, Re, Ee;

  R := E[1..3, 1..3];
  t := E[1..3, 4];
  te := <<round(evalf(t[1]*e^2))/e^2>,<round(evalf(t[2]*e^2))/e^2>,<round(evalf(t[3]*e^2))/e^2>>;
  q := R2q(R);
  qe := <<round(evalf(q[1]*e))/e>,<round(evalf(q[2]*e))/e>,<round(evalf(q[3]*e))/e>,<round(evalf(q[4]*e))/e>>;
  Re := q2R(qe);
  Ee := <Re|te>;
end proc:
