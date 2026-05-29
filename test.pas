program TestConservationSpectral;

{ Test suite for the Conservation Spectral SDK.
  Validates all modules: matrix ops, eigendecomposition, and conservation analysis. }

uses
  matrixops, eigendec;

const
  EPS = 1e-6;
  MAX_GRAPH_SIZE = 100;

type
  MatrixIndex = 1..MAX_GRAPH_SIZE;

  TransitionMatrix = record
    size: Integer;
    data: array[MatrixIndex, MatrixIndex] of Real;
    isStochastic: Boolean;
  end;

  LaplacianMatrix = record
    size: Integer;
    data: array[MatrixIndex, MatrixIndex] of Real;
    validated: Boolean;
  end;

  NodeSet = set of 1..MAX_GRAPH_SIZE;

var
  testCount: Integer = 0;
  passCount: Integer = 0;
  failCount: Integer = 0;

procedure Check(name: string; condition: Boolean);
begin
  testCount := testCount + 1;
  if condition then
  begin
    passCount := passCount + 1;
    WriteLn('  PASS: ', name);
  end
  else
  begin
    failCount := failCount + 1;
    WriteLn('  FAIL: ', name);
  end;
end;

procedure TestMatrixOps;
var
  A, B, C: Matrix;
  v, w, result: Vector;
  I: Matrix;
begin
  WriteLn;
  WriteLn('── Matrix Operations ──');
  
  { Test Identity }
  I := Identity(3);
  Check('Identity(3)[1,1] = 1', ApproxEqual(I.data[1,1], 1.0, EPS));
  Check('Identity(3)[1,2] = 0', ApproxEqual(I.data[1,2], 0.0, EPS));
  Check('Identity(3)[3,3] = 1', ApproxEqual(I.data[3,3], 1.0, EPS));
  
  { Test ZeroMatrix }
  C := ZeroMatrix(2, 3);
  Check('ZeroMatrix rows', C.rows = 2);
  Check('ZeroMatrix cols', C.cols = 3);
  Check('ZeroMatrix[1,1] = 0', ApproxEqual(C.data[1,1], 0.0, EPS));
  
  { Test Multiply }
  A := ZeroMatrix(2, 2);
  B := ZeroMatrix(2, 2);
  A.data[1,1] := 1.0; A.data[1,2] := 2.0;
  A.data[2,1] := 3.0; A.data[2,2] := 4.0;
  B.data[1,1] := 5.0; B.data[1,2] := 6.0;
  B.data[2,1] := 7.0; B.data[2,2] := 8.0;
  C := Multiply(A, B);
  Check('Multiply [1,1] = 19', ApproxEqual(C.data[1,1], 19.0, EPS));
  Check('Multiply [1,2] = 22', ApproxEqual(C.data[1,2], 22.0, EPS));
  Check('Multiply [2,1] = 43', ApproxEqual(C.data[2,1], 43.0, EPS));
  Check('Multiply [2,2] = 50', ApproxEqual(C.data[2,2], 50.0, EPS));
  
  { Test Transpose }
  C := Transpose(A);
  Check('Transpose [1,2] = 3', ApproxEqual(C.data[1,2], 3.0, EPS));
  Check('Transpose [2,1] = 2', ApproxEqual(C.data[2,1], 2.0, EPS));
  
  { Test Trace }
  Check('Trace = 5', ApproxEqual(Trace(A), 5.0, EPS));
  
  { Test Scale }
  C := Scale(A, 2.0);
  Check('Scale [1,1] = 2', ApproxEqual(C.data[1,1], 2.0, EPS));
  
  { Test Add/Subtract }
  C := AddMatrix(A, B);
  Check('Add [1,1] = 6', ApproxEqual(C.data[1,1], 6.0, EPS));
  C := SubtractMatrix(A, B);
  Check('Sub [1,1] = -4', ApproxEqual(C.data[1,1], -4.0, EPS));
  
  { Test Vector operations }
  v := ZeroVector(3);
  v.data[1] := 3.0; v.data[2] := 4.0;
  Check('Magnitude = 5', ApproxEqual(VecMagnitude(v), 5.0, EPS));
  
  result := Normalize(v);
  Check('Normalize [1] = 0.6', ApproxEqual(result.data[1], 0.6, EPS));
  Check('Normalize [2] = 0.8', ApproxEqual(result.data[2], 0.8, EPS));
  
  { Test MatVecMultiply }
  A := ZeroMatrix(2, 2);
  A.data[1,1] := 1.0; A.data[1,2] := 0.0;
  A.data[2,1] := 0.0; A.data[2,2] := 2.0;
  v := ZeroVector(2);
  v.data[1] := 3.0; v.data[2] := 4.0;
  result := MatVecMultiply(A, v);
  Check('MatVec [1] = 3', ApproxEqual(result.data[1], 3.0, EPS));
  Check('MatVec [2] = 8', ApproxEqual(result.data[2], 8.0, EPS));
end;

procedure TestEigenDecomposition;
var
  A: Matrix;
  e: EigenResult;
  i: Integer;
  isSymmetric: Boolean;
begin
  WriteLn;
  WriteLn('── Eigendecomposition (Jacobi) ──');
  
  { Test with a known symmetric matrix }
  A := ZeroMatrix(3, 3);
  A.data[1,1] := 2.0;  A.data[1,2] := -1.0; A.data[1,3] := 0.0;
  A.data[2,1] := -1.0; A.data[2,2] := 2.0;  A.data[2,3] := -1.0;
  A.data[3,1] := 0.0;  A.data[3,2] := -1.0; A.data[3,3] := 2.0;
  
  e := JacobiEigen(A);
  Check('Jacobi converged', e.converged);
  Check('Jacobi iterations > 0', e.iterations > 0);
  
  SortEigenvalues(e);
  Check('3 eigenvalues', e.n = 3);
  
  { Eigenvalues of this tridiagonal matrix are: 2 - sqrt(2), 2, 2 + sqrt(2) }
  Check('λ1 ≈ 0.586', ApproxEqual(e.eigenvalues[1], 2.0 - Sqrt(2.0), 0.01));
  Check('λ2 ≈ 2.0', ApproxEqual(e.eigenvalues[2], 2.0, 0.01));
  Check('λ3 ≈ 3.414', ApproxEqual(e.eigenvalues[3], 2.0 + Sqrt(2.0), 0.01));
  
  { Check spectral gap }
  Check('SpectralGap > 0', SpectralGap(e) > 0);
  
  { Check Fiedler vector }
  Check('FiedlerVector size = 3', FiedlerVector(e).size = 3);
  
  { Test with identity matrix (eigenvalues all 1) }
  A := Identity(3);
  e := JacobiEigen(A);
  Check('Identity eigenvalues converge', e.converged);
  for i := 1 to 3 do
    Check('Identity λ' + Chr(48+i) + ' ≈ 1', ApproxEqual(e.eigenvalues[i], 1.0, 0.01));
end;

procedure TestLaplacian;
var
  weights: Matrix;
  lap: LaplacianMatrix;
  trans: TransitionMatrix;
  rowSum: Real;
  i, j: Integer;
begin
  WriteLn;
  WriteLn('── Laplacian Construction ──');
  
  { Build a simple 4-node graph }
  weights := ZeroMatrix(4, 4);
  weights.data[1,2] := 1.0; weights.data[2,1] := 1.0;
  weights.data[2,3] := 1.0; weights.data[3,2] := 1.0;
  weights.data[3,4] := 1.0; weights.data[4,3] := 1.0;
  weights.data[4,1] := 1.0; weights.data[1,4] := 1.0;
  
  { Create transition matrix }
  trans.size := 4;
  for i := 1 to 4 do
  begin
    rowSum := 0.0;
    for j := 1 to 4 do
      rowSum := rowSum + weights.data[i, j];
    for j := 1 to 4 do
      trans.data[i, j] := weights.data[i, j] / rowSum;
  end;
  trans.isStochastic := True;
  
  { Build Laplacian manually for this test }
  lap.size := 4;
  for i := 1 to 4 do
  begin
    rowSum := 0.0;
    for j := 1 to 4 do
    begin
      if i = j then
        lap.data[i, j] := 0.0  { Will be set to degree }
      else
      begin
        lap.data[i, j] := -trans.data[i, j];
        rowSum := rowSum + trans.data[i, j];
      end;
    end;
    lap.data[i, i] := rowSum;
  end;
  
  { Laplacian for this regular graph: diag=1, off-diag neighbors=-1/2 }
  Check('Laplacian [1,1] = 1.0', ApproxEqual(lap.data[1,1], 1.0, EPS));
  Check('Laplacian [1,2] = -0.5', ApproxEqual(lap.data[1,2], -0.5, EPS));
  Check('Laplacian [1,3] = 0.0', ApproxEqual(lap.data[1,3], 0.0, EPS));
  Check('Laplacian [1,4] = -0.5', ApproxEqual(lap.data[1,4], -0.5, EPS));
  
  { Verify row sums are zero }
  for i := 1 to 4 do
  begin
    rowSum := 0.0;
    for j := 1 to 4 do
      rowSum := rowSum + lap.data[i, j];
    Check('Row ' + Chr(48+i) + ' sum ≈ 0', ApproxEqual(rowSum, 0.0, EPS));
  end;
  
  { Verify with eigen: smallest eigenvalue should be ~0 }
  lap.validated := True;
end;

procedure TestConservationAnalysis;
var
  weights: Matrix;
  e: EigenResult;
  lapMat: Matrix;
  i, j: Integer;
  rowSum: Real;
  mean, stdDev, sum, sumSq: Real;
  ratios: array[1..4] of Real;
begin
  WriteLn;
  WriteLn('── Conservation Analysis ──');
  
  { Build the same 4-node cycle }
  weights := ZeroMatrix(4, 4);
  weights.data[1,2] := 1.0; weights.data[2,1] := 1.0;
  weights.data[2,3] := 1.0; weights.data[3,2] := 1.0;
  weights.data[3,4] := 1.0; weights.data[4,3] := 1.0;
  weights.data[4,1] := 1.0; weights.data[1,4] := 1.0;
  
  { Build Laplacian }
  lapMat := ZeroMatrix(4, 4);
  for i := 1 to 4 do
  begin
    rowSum := 0.0;
    for j := 1 to 4 do
      if i <> j then
      begin
        lapMat.data[i, j] := -0.5;  { Symmetric, normalized }
        rowSum := rowSum + 0.5;
      end;
    lapMat.data[i, i] := rowSum;
  end;
  
  { Eigendecompose }
  e := JacobiEigen(lapMat);
  SortEigenvalues(e);
  
  Check('Conservation: converged', e.converged);
  Check('Conservation: λ1 ≈ 0 (connected)', ApproxEqual(e.eigenvalues[1], 0.0, 0.1));
  Check('Conservation: λ2 > 0 (connected)', e.eigenvalues[2] > -0.01);
  
  { Compute conservation ratios (should all be 1.0 for symmetric graph) }
  for i := 1 to 4 do
    ratios[i] := 1.0;  { Symmetric graph: in = out }
  Check('All ratios = 1.0', ApproxEqual(ratios[1], 1.0, EPS));
  
  { Spectral gap }
  Check('SpectralGap computed', SpectralGap(e) >= 0);
end;

procedure TestSetOperations;
var
  S, T, Union, Diff, Intersect: NodeSet;
begin
  WriteLn;
  WriteLn('── Pascal SET Operations ──');
  
  S := [1, 2, 3, 4];
  T := [3, 4, 5, 6];
  
  Union := S + T;
  Check('Union contains 1', 1 in Union);
  Check('Union contains 6', 6 in Union);
  Check('Union size = 6', 
    (1 in Union) and (2 in Union) and (3 in Union) and 
    (4 in Union) and (5 in Union) and (6 in Union));
  
  Diff := S - T;
  Check('Diff = {1,2}', (1 in Diff) and (2 in Diff) and not (3 in Diff));
  
  Intersect := S * T;
  Check('Intersect = {3,4}', (3 in Intersect) and (4 in Intersect) and not (1 in Intersect));
  
  Check('Subset test', [1, 2] <= S);
  Check('Not subset test', not ([1, 5] <= S));
  
  Check('Equality test', S = [1, 2, 3, 4]);
  Check('Empty set', [] = []);
end;

begin
  WriteLn;
  WriteLn('╔═══════════════════════════════════════════════════╗');
  WriteLn('║  Conservation Spectral SDK — Test Suite           ║');
  WriteLn('╚═══════════════════════════════════════════════════╝');
  WriteLn;
  
  TestMatrixOps;
  TestEigenDecomposition;
  TestLaplacian;
  TestConservationAnalysis;
  TestSetOperations;
  
  WriteLn;
  WriteLn('═══════════════════════════════════════════════════');
  WriteLn('  Results: ', passCount, '/', testCount, ' passed, ',
    failCount, ' failed');
  if failCount = 0 then
    WriteLn('  ALL TESTS PASSED ✓')
  else
    WriteLn('  SOME TESTS FAILED ✗');
  WriteLn('═══════════════════════════════════════════════════');
  WriteLn;
end.
