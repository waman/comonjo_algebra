# Development Notes (Japanese)
## 代数：algebra モジュール
coming soon...

## Polynomial 多項式：poly モジュール
この型は Scala ライブラリ Spire の [Polynomial](https://github.com/typelevel/spire/blob/main/core/src/main/scala/spire/math/Polynomial.scala) を概ね移植しています。

### 列挙型
多項式を表す構造体としては、最高次係数までの（0を含む）すべての係数を保持する **Dense** （密）と、0でない係数のみを次数をキーとして保持する **Sparse** （疎）が有用です。　Spire では、これらに加えて定数を表す Constant を加えて列挙型としていますが、comonjo_algebra では、定数の中でも 0 を別にして
- Zero （0次）
- Constant （1次）
- Dense （1次以上）
- Sparse （1次以上）

を `Polynomial` の列挙型にしています。 Zero を Constant と別にしているのは `num::traits::ConstZero` などとの兼ね合いもありますが、特に必須というわけではなく、実装の一選択です。

Dense では `Vec` を、Sparse では `BTreeMap` を用いて係数を保持しています（Spire では、Sparse を2つの配列で実装しています）。

### インスタンス生成
0 と定数を表す多項式は、`Polynomial` のファクトリ関数を用いて簡単にインスタンスを生成できます：
```rust
let zero = Polynomial::Zero();
let two = Polynomial::constant(2);
```
定数を生成する場合は、列挙型を直接生成するのではなく、型関連関数の `Polynomial::constant()` を呼び出す必要があります（定数が0だった場合に Zero を返すため）。

1次以上の多項式のオブジェクトを生成する場合、Dense か Sparse かを指定する必要があります。

#### マクロ
`dense![]` と `sparse![]` マクロを使って、高次の多項式を生成できます。　`dense![]` は昇順にすべての係数を列挙して引数に渡します。　また、`sparse![]` は項の次数 (`usize`) と係数のタプルを列挙して引数に渡します：
```rust
let p_d = dense![1, 2, 3];  // 1 + 2x + 3x²
let p_s = sparse![(0, 1), (2, 3)];  // 1 + 3x²
``` 
次数が0または1になる場合は、当然 `Polynomial::Zero()` もしくは `Polynomial::Constant()` が返されます。　`sparse![]` は、係数が0のタプルは無視します。

#### From, FromIterator
`Polynomial<C>` (where C: Semiring）は `From<Vec<C>>` と `From<BTreeMap<usize, C>>` を実装しているので、それぞれからインスタンス化できます：
```rust
// From<Vec<C>>
let v = vec![1, 2, 3];
let p_d = Polynomial::from(v);  // 1 + 2x + 3x²
assert_eq!(p_d, dense![1, 2, 3]);

// From<BTreeMap<usize, C>>
let mut m = BTreeMap::new();
m.insert(0, 1);
m.insert(2, 3);
let p_s = Polynomial::from(m);  // 1 + 3x²
assert_eq!(p_s, dense![1, 0, 3]);
```
同様に、`Polynomial<C>` は `FromIterator<C>` と `FromIterator<(usize, C)>` を実装しているので、以下のようにイテレーターの `collect()`でインスタンス化できます：
```rust
/// FromIterator<C>
let vec_iter = (1..8).into_iter();
let p: Polynomial<i64> = vec_iter.collect();
assert_eq!(p, dense![1, 2, 3, 4, 5, 6, 7]);

// FromIterator<(usize, C)>
let map_iter = (1..8).into_iter().enumerate()
                     .filter(|(i, c)|i % 2 == 0);
let p: Polynomial<i64> = map_iter.collect();
assert_eq!(p, sparse![(0, 1), (2, 3), (4, 5), (6, 7)]);
```

#### ファクトリ関数
#### parse
#### new_raw 関数 （パッケージ内）

### イテレーター
### 代数演算
### 数学関数