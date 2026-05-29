unit matrixops;

{ Matrix operations module for the Conservation Spectral SDK.
  Pascal's strong typing ensures every operation is type-safe. }

interface

const
  MAX_SIZE = 100;

type
  SizeRange = 1..MAX_SIZE;
  IndexRange = 1..MAX_SIZE;

  Matrix = record
    rows: Integer;
    cols: Integer;
    data: array[IndexRange, IndexRange] of Real;
  end;

  Vector = record
    size: Integer;
    data: array[IndexRange] of Real;
  end;

{ Create an identity matrix of given size }
function Identity(n: Integer): Matrix;

{ Create a zero matrix of given dimensions }
function ZeroMatrix(rows, cols: Integer): Matrix;

{ Create a zero vector of given size }
function ZeroVector(n: Integer): Vector;

{ Multiply two matrices: result = A * B }
function Multiply(const A, B: Matrix): Matrix;

{ Transpose a matrix }
function Transpose(const M: Matrix): Matrix;

{ Compute the trace of a square matrix }
function Trace(const M: Matrix): Real;

{ Multiply matrix by scalar }
function Scale(const M: Matrix; s: Real): Matrix;

{ Add two matrices }
function AddMatrix(const A, B: Matrix): Matrix;

{ Subtract: result = A - B }
function SubtractMatrix(const A, B: Matrix): Matrix;

{ Matrix-vector multiply: result = M * v }
function MatVecMultiply(const M: Matrix; const v: Vector): Vector;

{ Dot product of two vectors }
function DotProduct(const a, b: Vector): Vector;

{ Vector magnitude (L2 norm) }
function VecMagnitude(const v: Vector): Real;

{ Normalize a vector to unit length }
function Normalize(const v: Vector): Vector;

{ Deep copy a matrix }
function CopyMatrix(const M: Matrix): Matrix;

{ Print a matrix (for debugging) }
procedure PrintMatrix(const M: Matrix);

{ Print a vector }
procedure PrintVector(const v: Vector);

{ Check approximate equality of two reals }
function ApproxEqual(a, b: Real; eps: Real): Boolean;

implementation

function Identity(n: Integer): Matrix;
var
  i, j: Integer;
begin
  Identity.rows := n;
  Identity.cols := n;
  for i := 1 to MAX_SIZE do
    for j := 1 to MAX_SIZE do
      if (i = j) and (i <= n) then
        Identity.data[i, j] := 1.0
      else
        Identity.data[i, j] := 0.0;
end;

function ZeroMatrix(rows, cols: Integer): Matrix;
var
  i, j: Integer;
begin
  ZeroMatrix.rows := rows;
  ZeroMatrix.cols := cols;
  for i := 1 to MAX_SIZE do
    for j := 1 to MAX_SIZE do
      ZeroMatrix.data[i, j] := 0.0;
end;

function ZeroVector(n: Integer): Vector;
var
  i: Integer;
begin
  ZeroVector.size := n;
  for i := 1 to MAX_SIZE do
    ZeroVector.data[i] := 0.0;
end;

function Multiply(const A, B: Matrix): Matrix;
var
  i, j, k: Integer;
  sum: Real;
  result: Matrix;
begin
  result := ZeroMatrix(A.rows, B.cols);
  for i := 1 to A.rows do
    for j := 1 to B.cols do
    begin
      sum := 0.0;
      for k := 1 to A.cols do
        sum := sum + A.data[i, k] * B.data[k, j];
      result.data[i, j] := sum;
    end;
  Multiply := result;
end;

function Transpose(const M: Matrix): Matrix;
var
  i, j: Integer;
  result: Matrix;
begin
  result := ZeroMatrix(M.cols, M.rows);
  for i := 1 to M.rows do
    for j := 1 to M.cols do
      result.data[j, i] := M.data[i, j];
  Transpose := result;
end;

function Trace(const M: Matrix): Real;
var
  i: Integer;
  sum: Real;
begin
  sum := 0.0;
  for i := 1 to M.rows do
    if i <= M.cols then
      sum := sum + M.data[i, i];
  Trace := sum;
end;

function Scale(const M: Matrix; s: Real): Matrix;
var
  i, j: Integer;
  result: Matrix;
begin
  result := ZeroMatrix(M.rows, M.cols);
  for i := 1 to M.rows do
    for j := 1 to M.cols do
      result.data[i, j] := M.data[i, j] * s;
  Scale := result;
end;

function AddMatrix(const A, B: Matrix): Matrix;
var
  i, j: Integer;
  result: Matrix;
begin
  result := ZeroMatrix(A.rows, A.cols);
  for i := 1 to A.rows do
    for j := 1 to A.cols do
      result.data[i, j] := A.data[i, j] + B.data[i, j];
  AddMatrix := result;
end;

function SubtractMatrix(const A, B: Matrix): Matrix;
var
  i, j: Integer;
  result: Matrix;
begin
  result := ZeroMatrix(A.rows, A.cols);
  for i := 1 to A.rows do
    for j := 1 to A.cols do
      result.data[i, j] := A.data[i, j] - B.data[i, j];
  SubtractMatrix := result;
end;

function MatVecMultiply(const M: Matrix; const v: Vector): Vector;
var
  i, j: Integer;
  sum: Real;
  result: Vector;
begin
  result := ZeroVector(M.rows);
  for i := 1 to M.rows do
  begin
    sum := 0.0;
    for j := 1 to M.cols do
      sum := sum + M.data[i, j] * v.data[j];
    result.data[i] := sum;
  end;
  MatVecMultiply := result;
end;

function DotProduct(const a, b: Vector): Vector;
var
  i: Integer;
  sum: Real;
  result: Vector;
begin
  { Return a 1-element vector containing the dot product }
  result := ZeroVector(1);
  sum := 0.0;
  for i := 1 to a.size do
    sum := sum + a.data[i] * b.data[i];
  result.data[1] := sum;
  DotProduct := result;
end;

function VecMagnitude(const v: Vector): Real;
var
  i: Integer;
  sum: Real;
begin
  sum := 0.0;
  for i := 1 to v.size do
    sum := sum + v.data[i] * v.data[i];
  VecMagnitude := Sqrt(sum);
end;

function Normalize(const v: Vector): Vector;
var
  mag: Real;
  i: Integer;
  result: Vector;
begin
  mag := VecMagnitude(v);
  result := ZeroVector(v.size);
  if mag > 1e-10 then
    for i := 1 to v.size do
      result.data[i] := v.data[i] / mag;
  Normalize := result;
end;

function CopyMatrix(const M: Matrix): Matrix;
var
  i, j: Integer;
  result: Matrix;
begin
  result.rows := M.rows;
  result.cols := M.cols;
  for i := 1 to MAX_SIZE do
    for j := 1 to MAX_SIZE do
      result.data[i, j] := M.data[i, j];
  CopyMatrix := result;
end;

procedure PrintMatrix(const M: Matrix);
var
  i, j: Integer;
begin
  WriteLn('Matrix ', M.rows, 'x', M.cols, ':');
  for i := 1 to M.rows do
  begin
    for j := 1 to M.cols do
      Write(M.data[i, j]:10:4);
    WriteLn;
  end;
end;

procedure PrintVector(const v: Vector);
var
  i: Integer;
begin
  Write('Vector[', v.size, ']: ');
  for i := 1 to v.size do
    Write(v.data[i]:10:4);
  WriteLn;
end;

function ApproxEqual(a, b: Real; eps: Real): Boolean;
begin
  ApproxEqual := Abs(a - b) < eps;
end;

end.
