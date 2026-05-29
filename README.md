# Conservation Spectral SDK — Pascal Edition

**Type-safe structured spectral graph analysis, written in 1970s Pascal.**

This project demonstrates how Pascal's strong typing and structured programming constraints _improve_ algorithm design. The type system catches errors that would be silent bugs in C, and Pascal's unique SET type makes graph operations natural and elegant.

## What It Does

Given a weighted graph, the SDK:

1. **Builds a transition matrix** (row-stochastic, validated)
2. **Constructs the graph Laplacian** (type-safe, validated)
3. **Computes eigendecomposition** (Jacobi method, converged)
4. **Detects conservation anomalies** (nodes with unbalanced flow)
5. **Estimates the Cheeger constant** (spectral graph theory)

## Pascal's Unique Advantages

### 1. Type Safety
`TransitionMatrix`, `LaplacianMatrix`, and `ConservationReport` are **distinct types**. You cannot pass a Laplacian where a Transition is expected. The compiler catches this at compile time.

### 2. SET Type
Pascal has built-in set operations — no library needed:
```pascal
type
  NodeSet = set of 1..100;

function Neighbors(node: Integer; var graph: TransitionMatrix): NodeSet;
begin
  result := [];
  for j := 1 to graph.size do
    if graph.data[node, j] > 0 then
      result := result + [j];  { Set union! }
  Neighbors := result;
end;
```

### 3. Bounds Checking
Every array access is bounds-checked. No buffer overflows, no silent out-of-bounds reads.

### 4. Record Validation
The `LaplacianMatrix.validated` field is computed once and checked before use. The type system makes this a natural pattern.

### 5. Structured Thinking
Pascal forces you to declare types before using them. This upfront cost pays dividends in correctness and readability.

## Files

| File | Description |
|------|-------------|
| `conservation.pas` | Main program — graph types, Laplacian, conservation analysis |
| `matrixops.pas` | Matrix operations unit (multiply, transpose, etc.) |
| `eigendec.pas` | Jacobi eigendecomposition unit |
| `test.pas` | Test suite (55 tests) |

## Build & Run

Requires Free Pascal (`fpc`):

```bash
# Install (Ubuntu/Debian)
sudo apt install fp-compiler

# Compile units
fpc matrixops.pas
fpc eigendec.pas

# Compile and run main program
fpc conservation.pas
./conservation

# Compile and run tests
fpc test.pas
./test
```

## Sample Output

```
  Conservation Ratios (out/in weight ratio per node):
    Node 1:   0.7953
    Node 2:   0.7890
    Node 3:   0.8812
    Node 4:  14.4310  *** ANOMALY ***
    Node 5:   0.8846
    Node 6:   0.8766

  Spectral Gap:     0.209456
  Cheeger Constant (upper bound):     0.960370
```

Node 4 was deliberately given weak connections — the SDK detects it as a conservation anomaly.

## The Lesson

Pascal's constraints aren't limitations — they're _design tools_. By forcing you to think about types upfront, Pascal catches entire categories of bugs at compile time. The SET type, unique among mainstream languages, makes graph algorithms elegant. And the structured programming paradigm (no goto, nested functions, strong encapsulation) produces code that's readable decades later.

This is what happens when a language designed for _teaching_ meets a real algorithm: the code is better because the language won't let you be lazy.

## License

MIT
