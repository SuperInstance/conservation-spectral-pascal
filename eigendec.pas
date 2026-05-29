unit eigendec;

{ Jacobi eigendecomposition for symmetric matrices.
  Pascal's type system ensures we never mix up eigenvalues with eigenvectors. }

interface

uses
  matrixops;

const
  MAX_JACOBI_ITERATIONS = 1000;
  JACOBI_TOLERANCE = 1e-10;

type
  EigenResult = record
    n: Integer;
    eigenvalues: array[IndexRange] of Real;
    eigenvectors: Matrix;  { Columns are eigenvectors }
    iterations: Integer;
    converged: Boolean;
  end;

{ Perform Jacobi eigendecomposition on a symmetric matrix. }
function JacobiEigen(const A: Matrix): EigenResult;

{ Sort eigenvalues in ascending order (with corresponding eigenvectors) }
procedure SortEigenvalues(var e: EigenResult);

{ Compute the spectral gap }
function SpectralGap(const e: EigenResult): Real;

{ Compute Fiedler vector (eigenvector for 2nd smallest eigenvalue) }
function FiedlerVector(const e: EigenResult): Vector;

{ Return eigenvalue at given index }
function EigenvalueAt(const e: EigenResult; idx: Integer): Real;

implementation

function JacobiEigen(const A: Matrix): EigenResult;
var
  Work, V: Matrix;
  ii, jj, pp, qq, iter: Integer;
  maxOff, theta, cosVal, sinVal: Real;
  app, aqq, apq, aip, aiq, vip, viq: Real;
  eig: EigenResult;
  rowSum: Real;
begin
  Work := CopyMatrix(A);
  V := Identity(A.rows);

  eig.n := A.rows;
  eig.converged := False;
  eig.iterations := 0;

  for iter := 1 to MAX_JACOBI_ITERATIONS do
  begin
    eig.iterations := iter;

    { Find largest off-diagonal element }
    maxOff := 0.0;
    pp := 1;
    qq := 2;
    for ii := 1 to Work.rows do
      for jj := ii + 1 to Work.rows do
        if Abs(Work.data[ii, jj]) > maxOff then
        begin
          maxOff := Abs(Work.data[ii, jj]);
          pp := ii;
          qq := jj;
        end;

    if maxOff < JACOBI_TOLERANCE then
    begin
      eig.converged := True;
      Break;
    end;

    app := Work.data[pp, pp];
    aqq := Work.data[qq, qq];
    apq := Work.data[pp, qq];

    if Abs(app - aqq) < 1e-20 then
      theta := Pi / 4.0
    else
      theta := 0.5 * ArcTan(2.0 * apq / (app - aqq));

    cosVal := Cos(theta);
    sinVal := Sin(theta);

    { Apply rotation to Work }
    for ii := 1 to Work.rows do
    begin
      if (ii <> pp) and (ii <> qq) then
      begin
        aip := Work.data[ii, pp];
        aiq := Work.data[ii, qq];
        Work.data[ii, pp] := cosVal * aip + sinVal * aiq;
        Work.data[pp, ii] := Work.data[ii, pp];
        Work.data[ii, qq] := -sinVal * aip + cosVal * aiq;
        Work.data[qq, ii] := Work.data[ii, qq];
      end;
    end;

    Work.data[pp, pp] := cosVal * cosVal * app + 2 * sinVal * cosVal * apq + sinVal * sinVal * aqq;
    Work.data[qq, qq] := sinVal * sinVal * app - 2 * sinVal * cosVal * apq + cosVal * cosVal * aqq;
    Work.data[pp, qq] := 0.0;
    Work.data[qq, pp] := 0.0;

    { Update eigenvectors }
    for ii := 1 to V.rows do
    begin
      vip := V.data[ii, pp];
      viq := V.data[ii, qq];
      V.data[ii, pp] := cosVal * vip + sinVal * viq;
      V.data[ii, qq] := -sinVal * vip + cosVal * viq;
    end;
  end;

  for ii := 1 to Work.rows do
    eig.eigenvalues[ii] := Work.data[ii, ii];

  eig.eigenvectors := V;
  JacobiEigen := eig;
end;

procedure SortEigenvalues(var e: EigenResult);
var
  ii, jj, minIdx: Integer;
  tmpVal: Real;
begin
  for ii := 1 to e.n - 1 do
  begin
    minIdx := ii;
    for jj := ii + 1 to e.n do
      if e.eigenvalues[jj] < e.eigenvalues[minIdx] then
        minIdx := jj;

    if minIdx <> ii then
    begin
      tmpVal := e.eigenvalues[ii];
      e.eigenvalues[ii] := e.eigenvalues[minIdx];
      e.eigenvalues[minIdx] := tmpVal;

      for jj := 1 to e.n do
      begin
        tmpVal := e.eigenvectors.data[jj, ii];
        e.eigenvectors.data[jj, ii] := e.eigenvectors.data[jj, minIdx];
        e.eigenvectors.data[jj, minIdx] := tmpVal;
      end;
    end;
  end;
end;

function SpectralGap(const e: EigenResult): Real;
var
  ii: Integer;
  prev: Real;
begin
  prev := 0.0;
  for ii := 1 to e.n do
  begin
    if e.eigenvalues[ii] > JACOBI_TOLERANCE then
    begin
      if prev > JACOBI_TOLERANCE then
      begin
        SpectralGap := e.eigenvalues[ii] - prev;
        Exit;
      end;
      prev := e.eigenvalues[ii];
    end;
  end;
  SpectralGap := 0.0;
end;

function FiedlerVector(const e: EigenResult): Vector;
var
  ii: Integer;
  vec: Vector;
begin
  vec := ZeroVector(e.n);
  for ii := 1 to e.n do
    vec.data[ii] := e.eigenvectors.data[ii, 2];
  FiedlerVector := vec;
end;

function EigenvalueAt(const e: EigenResult; idx: Integer): Real;
begin
  if (idx >= 1) and (idx <= e.n) then
    EigenvalueAt := e.eigenvalues[idx]
  else
    EigenvalueAt := 0.0;
end;

end.
