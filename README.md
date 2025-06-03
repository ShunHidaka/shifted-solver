# Krylov subspace methods for shifted linear systems
Solver library for the standard shifted linear systems:
```math
(A + \sigma^{(m)} I) \textbf{x}^{(m)} = \textbf{b},\qquad (m=1,\dots,M),
```
* shifted CG method  
  for real symmetric / complex Hermitian **positive definite** matrix  + real shifts
* shifted MINRES method  
  for real symmetric / complex Hermitian matrix + real / complex shifts

## 更新履歴
2025/06/03: shifted CG法にseed switchingを実装

## TODO
* ドキュメントを作成する
* もう少し大きい行列での動作の確認
* sminresのGivens回転で自作関数を使用しているので他にないか探す
* 高速化・最適化