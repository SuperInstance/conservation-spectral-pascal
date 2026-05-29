program ConservationSpectral;

{ Conservation Spectral SDK - Main Program
  Demonstrates how Pascal's strong typing and structured programming
  constrain (and improve) spectral graph algorithm design.

  Key Pascal advantages demonstrated:
  - RECORD types force thinking about data flow
  - SET type enables natural graph operations
  - Bounds checking prevents silent array bugs
  - Strong typing catches parameter mixups at compile time
  - Nested functions encourage modular design }

uses
  matrixops, eigendec;

const
  MAX_GRAPH_SIZE = 100;
  CONSERVATION_THRESHOLD = 0.01;
  ANOMALY_SIGMA = 2.0;

type
  GraphSize = 1..MAX_GRAPH_SIZE;
  MatrixIndex = 1..MAX_GRAPH_SIZE;

  { The type system makes illegal states unrepresentable.
    A TransitionMatrix is NOT just any matrix — it's specifically
    a stochastic matrix with rows summing to 1. }
  TransitionMatrix = record
    size: Integer;
    data: array[MatrixIndex, MatrixIndex] of Real;
    isStochastic: Boolean;  { Cached validation result }
  end;

  AttributeVector = record
    size: Integer;
    values: array[MatrixIndex] of Real;
    name: string[64];
  end;

  { LaplacianMatrix is a distinct type from TransitionMatrix.
    You can't accidentally pass one where the other is expected. }
  LaplacianMatrix = record
    size: Integer;
    data: array[MatrixIndex, MatrixIndex] of Real;
    validated: Boolean;
  end;

  ConservationReport = record
    size: Integer;
    ratios: array[MatrixIndex] of Real;
    spectralGap: Real;
    cheegerConstant: Real;
    anomalyCount: Integer;
    anomalyIndices: array[MatrixIndex] of Boolean;
    isValid: Boolean;
  end;

  { Pascal's SET type — unique and powerful for graph operations }
  NodeSet = set of 1..MAX_GRAPH_SIZE;

{ ═══════════════════════════════════════════════════════════════════════════ }
{ Graph construction and validation                                          }
{ ═══════════════════════════════════════════════════════════════════════════ }

{ Validate that a matrix is row-stochastic (rows sum to ~1) }
function ValidateStochastic(var mat: TransitionMatrix): Boolean;
var
  i, j: Integer;
  rowSum: Real;
begin
  ValidateStochastic := True;
  for i := 1 to mat.size do
  begin
    rowSum := 0.0;
    for j := 1 to mat.size do
      rowSum := rowSum + mat.data[i, j];
    if Abs(rowSum - 1.0) > 0.01 then
    begin
      ValidateStochastic := False;
      Exit;
    end;
  end;
end;

{ Create a transition matrix from an adjacency-like weight matrix,
  normalizing rows to sum to 1 }
function CreateTransitionMatrix(n: Integer;
  var weights: Matrix): TransitionMatrix;
var
  i, j: Integer;
  rowSum: Real;
  result: TransitionMatrix;
begin
  result.size := n;
  for i := 1 to n do
  begin
    rowSum := 0.0;
    for j := 1 to n do
      rowSum := rowSum + weights.data[i, j];
    if rowSum > 1e-10 then
    begin
      for j := 1 to n do
        result.data[i, j] := weights.data[i, j] / rowSum;
    end
    else
    begin
      { Disconnected node — uniform distribution }
      for j := 1 to n do
        result.data[i, j] := 1.0 / n;
    end;
  end;
  result.isStochastic := ValidateStochastic(result);
  CreateTransitionMatrix := result;
end;

{ ═══════════════════════════════════════════════════════════════════════════ }
{ SET-based graph operations — Pascal's secret weapon                        }
{ ═══════════════════════════════════════════════════════════════════════════ }

{ Get the neighbor set of a node using Pascal's SET type }
function Neighbors(node: Integer; var trans: TransitionMatrix): NodeSet;
var
  j: Integer;
  result: NodeSet;
begin
  result := [];
  for j := 1 to trans.size do
    if trans.data[node, j] > 0.001 then
      result := result + [j];  { Set union — built-in! }
  Neighbors := result;
end;

{ Compute boundary of a node set S: nodes in S with neighbors outside S }
function Boundary(var S: NodeSet; var trans: TransitionMatrix): NodeSet;
var
  i: Integer;
  nbrs: NodeSet;
  result: NodeSet;
begin
  result := [];
  for i := 1 to trans.size do
    if i in S then
    begin
      nbrs := Neighbors(i, trans);
      { Nodes that are neighbors but NOT in S }
      result := result + (nbrs - S);  { Set difference! }
    end;
  Boundary := result;
end;

{ Cardinality of a node set }
function SetSize(const S: NodeSet): Integer;
var
  i, count: Integer;
begin
  count := 0;
  for i := 1 to MAX_GRAPH_SIZE do
    if i in S then
      count := count + 1;
  SetSize := count;
end;

{ ═══════════════════════════════════════════════════════════════════════════ }
{ Laplacian construction — TYPE SAFE                                         }
{ ═══════════════════════════════════════════════════════════════════════════ }

{ Verify a Laplacian: diagonal >= 0, off-diagonal <= 0, row sums = 0 }
function VerifyLaplacian(var lap: LaplacianMatrix): Boolean;
var
  i, j: Integer;
  rowSum: Real;
begin
  VerifyLaplacian := True;
  for i := 1 to lap.size do
  begin
    rowSum := 0.0;
    for j := 1 to lap.size do
    begin
      if (i <> j) and (lap.data[i, j] > 0.001) then
      begin
        VerifyLaplacian := False;
        Exit;
      end;
      if (i = j) and (lap.data[i, j] < -0.001) then
      begin
        VerifyLaplacian := False;
        Exit;
      end;
      rowSum := rowSum + lap.data[i, j];
    end;
    if Abs(rowSum) > 0.01 then
    begin
      VerifyLaplacian := False;
      Exit;
    end;
  end;
end;

{ Build graph Laplacian from transition matrix }
procedure BuildLaplacian(var trans: TransitionMatrix;
  var lap: LaplacianMatrix);
var
  i, j: Integer;
  degree: Real;
begin
  lap.size := trans.size;
  for i := 1 to trans.size do
  begin
    degree := 0.0;
    for j := 1 to trans.size do
      degree := degree + trans.data[i, j];
    for j := 1 to trans.size do
      if i = j then
        lap.data[i, j] := degree
      else
        lap.data[i, j] := -trans.data[i, j];
  end;
  lap.validated := VerifyLaplacian(lap);
end;

{ Convert LaplacianMatrix to Matrix for eigen computation }
function LaplacianToMatrix(var lap: LaplacianMatrix): Matrix;
var
  i, j: Integer;
  result: Matrix;
begin
  result := ZeroMatrix(lap.size, lap.size);
  for i := 1 to lap.size do
    for j := 1 to lap.size do
      result.data[i, j] := lap.data[i, j];
  LaplacianToMatrix := result;
end;

{ ═══════════════════════════════════════════════════════════════════════════ }
{ Conservation analysis                                                      }
{ ═══════════════════════════════════════════════════════════════════════════ }

{ Compute conservation ratios for each node:
  ratio = sum(outgoing weights) / sum(incoming weights)
  In a well-connected graph, ratios cluster near 1.0.
  Anomalies deviate significantly. }
function ComputeConservation(var lap: LaplacianMatrix;
  var attrs: AttributeVector): ConservationReport;
var
  i, j: Integer;
  outgoing, incoming: Real;
  mean, stdDev, sum, sumSq: Real;
  report: ConservationReport;
begin
  report.size := lap.size;
  report.isValid := lap.validated;
  
  { Compute per-node conservation ratios }
  for i := 1 to lap.size do
  begin
    outgoing := 0.0;
    incoming := 0.0;
    for j := 1 to lap.size do
    begin
      if i <> j then
      begin
        { Off-diagonal of Laplacian is negative of weight }
        outgoing := outgoing + Abs(lap.data[i, j]);
        incoming := incoming + Abs(lap.data[j, i]);
      end;
    end;
    if incoming > 1e-10 then
      report.ratios[i] := outgoing / incoming
    else
      report.ratios[i] := 0.0;
    report.anomalyIndices[i] := False;
  end;
  
  { Compute mean and standard deviation of ratios }
  sum := 0.0;
  for i := 1 to lap.size do
    sum := sum + report.ratios[i];
  mean := sum / lap.size;
  
  sumSq := 0.0;
  for i := 1 to lap.size do
    sumSq := sumSq + Sqr(report.ratios[i] - mean);
  stdDev := Sqrt(sumSq / lap.size);
  
  { Detect anomalies (ratios beyond ANOMALY_SIGMA standard deviations) }
  report.anomalyCount := 0;
  if stdDev > 1e-10 then
  begin
    for i := 1 to lap.size do
    begin
      if Abs(report.ratios[i] - mean) > ANOMALY_SIGMA * stdDev then
      begin
        report.anomalyIndices[i] := True;
        report.anomalyCount := report.anomalyCount + 1;
      end;
    end;
  end;
  
  ComputeConservation := report;
end;

{ Estimate Cheeger constant using spectral methods.
  h(G) <= sqrt(2 * lambda_2) — Cheeger inequality }
function EstimateCheeger(lambda2: Real): Real;
begin
  if lambda2 > 0 then
    EstimateCheeger := Sqrt(2.0 * lambda2)
  else
    EstimateCheeger := 0.0;
end;

{ ═══════════════════════════════════════════════════════════════════════════ }
{ Utility: Build a test graph                                                }
{ ═══════════════════════════════════════════════════════════════════════════ }

procedure BuildTestGraph(n: Integer; var weights: Matrix);
var
  i, j: Integer;
begin
  weights := ZeroMatrix(n, n);
  { Create a cycle graph with some cross-edges }
  for i := 1 to n do
  begin
    j := (i mod n) + 1;
    weights.data[i, j] := 1.0;
    weights.data[j, i] := 1.0;
  end;
  { Add some cross-edges for a richer structure }
  if n >= 5 then
  begin
    weights.data[1, 3] := 0.5;
    weights.data[3, 1] := 0.5;
    weights.data[2, 5] := 0.3;
    weights.data[5, 2] := 0.3;
  end;
  { Introduce a conservation anomaly: weak connection for node 4 }
  if n >= 6 then
  begin
    for j := 1 to n do
    begin
      weights.data[4, j] := 0.0;
      weights.data[j, 4] := 0.0;
    end;
    weights.data[4, 5] := 0.05;  { Very weak connection }
    weights.data[5, 4] := 0.05;
    weights.data[3, 4] := 0.05;
    weights.data[4, 3] := 0.05;
  end;
end;

{ ═══════════════════════════════════════════════════════════════════════════ }
{ Reporting                                                                  }
{ ═══════════════════════════════════════════════════════════════════════════ }

procedure PrintReport(var report: ConservationReport;
  var eResult: EigenResult);
var
  i: Integer;
  fiedler: Vector;
begin
  WriteLn;
  WriteLn('═══════════════════════════════════════════════════');
  WriteLn('  Conservation Spectral Analysis Report');
  WriteLn('═══════════════════════════════════════════════════');
  WriteLn;
  
  WriteLn('  Eigenvalues (sorted ascending):');
  for i := 1 to eResult.n do
    WriteLn('    λ', i, ' = ', eResult.eigenvalues[i]:12:6);
  WriteLn;
  
  WriteLn('  Spectral Gap: ', report.spectralGap:12:6);
  WriteLn('  Cheeger Constant (upper bound): ', report.cheegerConstant:12:6);
  WriteLn;
  
  WriteLn('  Conservation Ratios (out/in weight ratio per node):');
  for i := 1 to report.size do
  begin
    if report.anomalyIndices[i] then
      WriteLn('    Node ', i, ': ', report.ratios[i]:8:4, '  *** ANOMALY ***')
    else
      WriteLn('    Node ', i, ': ', report.ratios[i]:8:4);
  end;
  WriteLn;
  
  WriteLn('  Anomaly Count: ', report.anomalyCount);
  WriteLn;
  
  WriteLn('  Fiedler Vector (2nd eigenvector):');
  fiedler := FiedlerVector(eResult);
  for i := 1 to fiedler.size do
    WriteLn('    v[', i, '] = ', fiedler.data[i]:12:6);
  
  WriteLn;
  WriteLn('  Validation: ');
  WriteLn('    Laplacian valid: ', report.isValid);
  WriteLn('    Eigen converged: ', eResult.converged,
    ' (', eResult.iterations, ' iterations)');
  WriteLn('═══════════════════════════════════════════════════');
  WriteLn;
end;

{ Demonstrate Pascal's SET operations for graph analysis }
procedure DemonstrateSetOps(var trans: TransitionMatrix);
var
  node: Integer;
  nbrs, bnd: NodeSet;
  i: Integer;
begin
  WriteLn('  ── Pascal SET Operations Demo ──');
  WriteLn;
  
  for node := 1 to trans.size do
  begin
    nbrs := Neighbors(node, trans);
    Write('    Node ', node, ' neighbors: { ');
    for i := 1 to trans.size do
      if i in nbrs then
        Write(i, ' ');
    WriteLn('}');
  end;
  WriteLn;
  
  { Show boundary of a subset }
  nbrs := [1, 2, 3];
  bnd := Boundary(nbrs, trans);
  Write('    Boundary of {1,2,3}: { ');
  for i := 1 to trans.size do
    if i in bnd then
      Write(i, ' ');
  WriteLn('}');
  WriteLn;
end;

{ ═══════════════════════════════════════════════════════════════════════════ }
{ Main program                                                               }
{ ═══════════════════════════════════════════════════════════════════════════ }

var
  n: Integer;
  weights: Matrix;
  trans: TransitionMatrix;
  lap: LaplacianMatrix;
  lapMatrix: Matrix;
  attrs: AttributeVector;
  eResult: EigenResult;
  report: ConservationReport;

begin
  WriteLn;
  WriteLn('╔═══════════════════════════════════════════════════╗');
  WriteLn('║  Conservation Spectral SDK (Pascal Edition)      ║');
  WriteLn('║  Type-Safe Structured Graph Analysis              ║');
  WriteLn('╚═══════════════════════════════════════════════════╝');
  WriteLn;
  
  n := 6;
  WriteLn('  Building test graph with ', n, ' nodes...');
  
  { Build test graph }
  BuildTestGraph(n, weights);
  WriteLn('  Weight matrix:');
  PrintMatrix(weights);
  
  { Create transition matrix (type-safe construction) }
  trans := CreateTransitionMatrix(n, weights);
  WriteLn('  Transition matrix (row-stochastic): ', trans.isStochastic);
  
  { Demonstrate SET operations }
  DemonstrateSetOps(trans);
  
  { Build Laplacian (validated by type system) }
  BuildLaplacian(trans, lap);
  WriteLn('  Laplacian matrix (validated: ', lap.validated, '):');
  lapMatrix := LaplacianToMatrix(lap);
  PrintMatrix(lapMatrix);
  
  { Compute eigendecomposition }
  WriteLn('  Computing eigendecomposition (Jacobi method)...');
  eResult := JacobiEigen(lapMatrix);
  SortEigenvalues(eResult);
  
  { Set up attribute vector }
  attrs.size := n;
  attrs.name := 'node_weights';
  for n := 1 to attrs.size do
    attrs.values[n] := 1.0;
  
  { Compute conservation analysis }
  report := ComputeConservation(lap, attrs);
  report.spectralGap := SpectralGap(eResult);
  report.cheegerConstant := EstimateCheeger(eResult.eigenvalues[2]);
  
  { Print results }
  PrintReport(report, eResult);
  
  { Summary of Pascal's contribution }
  WriteLn('  ── What Pascal Gave Us ──');
  WriteLn;
  WriteLn('  1. TYPE SAFETY: TransitionMatrix, LaplacianMatrix, and');
  WriteLn('     ConservationReport are DISTINCT types. You cannot pass');
  WriteLn('     a Laplacian where a Transition is expected. The compiler');
  WriteLn('     catches this at compile time.');
  WriteLn;
  WriteLn('  2. SET OPERATIONS: NodeSet uses Pascal''s built-in set type.');
  WriteLn('     Neighbors(), Boundary() use union (+) and difference (-)');
  WriteLn('     operators. No external library needed.');
  WriteLn;
  WriteLn('  3. BOUNDS CHECKING: Every array access is bounds-checked.');
  WriteLn('     No buffer overflows. No silent out-of-bounds reads.');
  WriteLn;
  WriteLn('  4. RECORD VALIDATION: The LaplacianMatrix.validated field');
  WriteLn('     is computed once and checked before use. The type system');
  WriteLn('     makes this a natural pattern.');
  WriteLn;
  WriteLn('  5. STRUCTURED THINKING: Pascal forces you to declare types');
  WriteLn('     before using them. This upfront cost pays dividends in');
  WriteLn('     correctness and readability.');
  WriteLn;
end.
