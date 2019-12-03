with(ListTools):
with(LinearAlgebra):
with(PolynomialTools):
with(combinat, choose):
with(Groebner):
with(MatrixPolynomialAlgebra):
interface(rtablesize=24):
interface(warnlevel=0):
Digits:=30:
RatResolution:=10000:

kernelopts(numcpus=1):

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
  <c,s>;
end proc:

E2ER:=proc(E,e)
  local R, t;
  local te, q, qe, Re, Ee;

  R := E[1..3, 1..3];
  t := E[1..3, 4];
  te := <<round(evalf(t[1]*e))/e>,<round(evalf(t[2]*e))/e>,<round(evalf(t[3]*e))/e>>;
  q := R2q(R);
  qe := <<round(evalf(q[1]*e))/e>,<round(evalf(q[2]*e))/e>,<round(evalf(q[3]*e))/e>,<round(evalf(q[4]*e))/e>>;
  Re := q2R(qe);
  Ee := <Re|te>;
end proc:

# mechanism description
MechParams := {
a1=0,    a2=0, a3=0,     a4=0,    a5=0,     a6=0, a7=0,
d1=340,   d2=0,  d3=400,     d4=0,  d5=400,     d6=0, d7=126,
A1=-Pi/2, A2=Pi/2,   A3=-Pi/2, A4=Pi/2, A5=-Pi/2, A6=Pi/2, A7=0};

Mechanism := {op(MechParams),
              lambda1=CS2CSR(subs(MechParams,A1),RatResolution)[1],mu1=CS2CSR(subs(MechParams,A1),RatResolution)[2],
              lambda2=CS2CSR(subs(MechParams,A2),RatResolution)[1],mu2=CS2CSR(subs(MechParams,A2),RatResolution)[2],
              lambda3=CS2CSR(subs(MechParams,A3),RatResolution)[1],mu3=CS2CSR(subs(MechParams,A3),RatResolution)[2],
              lambda4=CS2CSR(subs(MechParams,A4),RatResolution)[1],mu4=CS2CSR(subs(MechParams,A4),RatResolution)[2],
              lambda5=CS2CSR(subs(MechParams,A5),RatResolution)[1],mu5=CS2CSR(subs(MechParams,A5),RatResolution)[2],
              lambda6=CS2CSR(subs(MechParams,A6),RatResolution)[1],mu6=CS2CSR(subs(MechParams,A6),RatResolution)[2],
              lambda7=CS2CSR(subs(MechParams,A7),RatResolution)[1],mu7=CS2CSR(subs(MechParams,A7),RatResolution)[2]};

# import end effector pose
Mh := Matrix(4,4,[[lx,mx,nx,rx],[ly,my,ny,ry],[lz,mz,nz,rz],[0,0,0,1]]):
MhV := ImportMatrix(MhFile)[2][1..3, 1..4];
MhR := E2ER(MhV,RatResolution);
MhRV := map(x->x[1]=x[2],convert(<convert(Mh[1..3,1..4],Vector)|convert(MhR,Vector)>,listlist));

# equations
lrE:=M2L((dhT(3).dhT(4).dhT(5)- dhInv(dhT(2)).dhInv(dhT(1)).<<MhR>,<0|0|0|1>>.dhInv(dhT(7)).dhInv(dhT(6)))[1..3,1..4]):
eq:=simplify(subs([op(convert(evalf(Mechanism),rational)),op(MhRV)],M2L(lrE)));

eqX:=numer~(normal~(([op(eq),s1^2+c1^2-1,s2^2+c2^2-1,s3^2+c3^2-1,s4^2+c4^2-1,s5^2+c5^2-1,s6^2+c6^2-1,s7^2+c7^2-1]))):

infolevel[GroebnerBasis]:=0:
B:=Basis(eqX, tdeg(op(indets(eqX)))):

e := eliminate(B, [c4,s4]):
vrEl := e[1];
eqXEl := e[2]:
vrNonel:=[c1,s1,c2,s2,c3,s3,c5,s5,c6,s6,c7,s7]:
BEl:=Basis(eqXEl,tdeg(op(indets(eqXEl)))):

idxDegTwo := ListTools[SearchAll](2, map(degree, BEl)):
eqQ := BEl[[idxDegTwo]]:

BElNew := Basis(eqQ, tdeg(op(indets(eqXEl)))):

with(CodeGeneration):
if evalb(BEl = BElNew) then

  for i from 1 to nops(eqQ) do
    eqQ[i]:= evalf(eqQ[i]/LeadingCoefficient(eqQ[i], tdeg(op(vrNonel)))):
  end do:
  eqMatlab := Matlab(eqQ, output=string):
  fd := fopen(eqFile, WRITE, TEXT):
  FileTools[Text][WriteString](fd, eqMatlab):
  FileTools[Text][WriteString](fd, Matlab(evalf(op(vrEl[1])[2]), output=string)):
  FileTools[Text][WriteString](fd, Matlab(simplify(evalf(op(vrEl[2])[2])), output=string)):
  fclose(fd):
end if;
