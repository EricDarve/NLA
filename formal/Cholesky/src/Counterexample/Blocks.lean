import src.Norms

/-!
# Block matrices indexed by `Fin (a + b)`

`blk A B C D` is the block matrix `[[A, B], [C, D]]`, with the first block indexed by
`Fin.castAdd` and the second by `Fin.natAdd`. This file collects the block formulas for
entries, products, transposes, matrix-vector products, and Euclidean norms.
-/
namespace Cholesky
open Matrix

variable {a b c d e f : ℕ}

/-- The block matrix `[[A, B], [C, D]]`. -/
def blk (A : Matrix (Fin a) (Fin c) ℝ) (B : Matrix (Fin a) (Fin d) ℝ)
    (C : Matrix (Fin b) (Fin c) ℝ) (D : Matrix (Fin b) (Fin d) ℝ) :
    Matrix (Fin (a + b)) (Fin (c + d)) ℝ :=
  Matrix.of fun i j =>
    Fin.addCases (motive := fun _ => ℝ) (fun i1 => Fin.addCases (motive := fun _ => ℝ) (A i1) (B i1) j)
      (fun i2 => Fin.addCases (motive := fun _ => ℝ) (C i2) (D i2) j) i

section
variable (A : Matrix (Fin a) (Fin c) ℝ) (B : Matrix (Fin a) (Fin d) ℝ)
  (C : Matrix (Fin b) (Fin c) ℝ) (D : Matrix (Fin b) (Fin d) ℝ)

@[simp] lemma blk_ll (i : Fin a) (j : Fin c) :
    blk A B C D (Fin.castAdd b i) (Fin.castAdd d j) = A i j := by simp [blk]
@[simp] lemma blk_lr (i : Fin a) (j : Fin d) :
    blk A B C D (Fin.castAdd b i) (Fin.natAdd c j) = B i j := by simp [blk]
@[simp] lemma blk_rl (i : Fin b) (j : Fin c) :
    blk A B C D (Fin.natAdd a i) (Fin.castAdd d j) = C i j := by simp [blk]
@[simp] lemma blk_rr (i : Fin b) (j : Fin d) :
    blk A B C D (Fin.natAdd a i) (Fin.natAdd c j) = D i j := by simp [blk]

lemma blk_transpose : (blk A B C D)ᵀ = blk Aᵀ Cᵀ Bᵀ Dᵀ := by
  ext i j
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> simp

/-- The first and second parts of a vector on `Fin (a + b)`. -/
def vL (v : Fin (a + b) → ℝ) : Fin a → ℝ := fun i => v (Fin.castAdd b i)
def vR (v : Fin (a + b) → ℝ) : Fin b → ℝ := fun i => v (Fin.natAdd a i)

lemma append_vL_vR (v : Fin (a + b) → ℝ) : Fin.append (vL v) (vR v) = v :=
  Fin.append_castAdd_natAdd

@[simp] lemma vL_append (x : Fin a → ℝ) (y : Fin b → ℝ) : vL (Fin.append x y) = x := by
  ext i; simp [vL]
@[simp] lemma vR_append (x : Fin a → ℝ) (y : Fin b → ℝ) : vR (Fin.append x y) = y := by
  ext i; simp [vR]

lemma blk_mulVec (v : Fin (c + d) → ℝ) :
    blk A B C D *ᵥ v = Fin.append (A *ᵥ vL v + B *ᵥ vR v) (C *ᵥ vL v + D *ᵥ vR v) := by
  ext i
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_add, vL, vR]

lemma blk_mulVec_append (x : Fin c → ℝ) (y : Fin d → ℝ) :
    blk A B C D *ᵥ Fin.append x y = Fin.append (A *ᵥ x + B *ᵥ y) (C *ᵥ x + D *ᵥ y) := by
  rw [blk_mulVec]; simp

end

lemma append_dotProduct (x z : Fin a → ℝ) (y w : Fin b → ℝ) :
    Fin.append x y ⬝ᵥ Fin.append z w = x ⬝ᵥ z + y ⬝ᵥ w := by
  simp [dotProduct, Fin.sum_univ_add]

lemma vnorm_append_sq (x : Fin a → ℝ) (y : Fin b → ℝ) :
    vnorm (Fin.append x y) ^ 2 = vnorm x ^ 2 + vnorm y ^ 2 := by
  rw [vnorm_sq, vnorm_sq, vnorm_sq, append_dotProduct]

lemma vnorm_sq_split (v : Fin (a + b) → ℝ) : vnorm v ^ 2 = vnorm (vL v) ^ 2 + vnorm (vR v) ^ 2 := by
  conv_lhs => rw [← append_vL_vR v]
  exact vnorm_append_sq _ _

lemma append_add (x x' : Fin a → ℝ) (y y' : Fin b → ℝ) :
    Fin.append x y + Fin.append x' y' = Fin.append (x + x') (y + y') := by
  ext i; refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;> simp

lemma append_smul (r : ℝ) (x : Fin a → ℝ) (y : Fin b → ℝ) :
    r • Fin.append x y = Fin.append (r • x) (r • y) := by
  ext i; refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;> simp

lemma append_zero : Fin.append (0 : Fin a → ℝ) (0 : Fin b → ℝ) = 0 := by
  ext i; refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;> simp

lemma append_eq_zero_iff (x : Fin a → ℝ) (y : Fin b → ℝ) :
    Fin.append x y = 0 ↔ x = 0 ∧ y = 0 := by
  constructor
  · intro h
    constructor
    · ext i; have := congrFun h (Fin.castAdd b i); simpa using this
    · ext i; have := congrFun h (Fin.natAdd a i); simpa using this
  · rintro ⟨rfl, rfl⟩; exact append_zero

/-- Block multiplication. -/
lemma blk_mul (A : Matrix (Fin a) (Fin c) ℝ) (B : Matrix (Fin a) (Fin d) ℝ)
    (C : Matrix (Fin b) (Fin c) ℝ) (D : Matrix (Fin b) (Fin d) ℝ)
    (A' : Matrix (Fin c) (Fin e) ℝ) (B' : Matrix (Fin c) (Fin f) ℝ)
    (C' : Matrix (Fin d) (Fin e) ℝ) (D' : Matrix (Fin d) (Fin f) ℝ) :
    blk A B C D * blk A' B' C' D' =
      blk (A * A' + B * C') (A * B' + B * D') (C * A' + D * C') (C * B' + D * D') := by
  ext i j
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_add]

/-- All-ones matrix. -/
def ones (a b : ℕ) : Matrix (Fin a) (Fin b) ℝ := Matrix.of fun _ _ => 1

lemma ones_mulVec (v : Fin b → ℝ) : ones a b *ᵥ v = fun _ => ∑ j, v j := by
  ext i; simp [ones, Matrix.mulVec, dotProduct]

lemma dot_ones_mulVec (w : Fin a → ℝ) (v : Fin b → ℝ) :
    w ⬝ᵥ (ones a b *ᵥ v) = (∑ i, w i) * ∑ j, v j := by
  rw [ones_mulVec]; simp [dotProduct, Finset.sum_mul]

lemma ones_eq_blk : ones (a + b) (c + d) = blk (ones a c) (ones a d) (ones b c) (ones b d) := by
  ext i j
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> simp [ones]

lemma blk_add (A A' : Matrix (Fin a) (Fin c) ℝ) (B B' : Matrix (Fin a) (Fin d) ℝ)
    (C C' : Matrix (Fin b) (Fin c) ℝ) (D D' : Matrix (Fin b) (Fin d) ℝ) :
    blk A B C D + blk A' B' C' D' = blk (A + A') (B + B') (C + C') (D + D') := by
  ext i j
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> simp

lemma smul_blk (r : ℝ) (A : Matrix (Fin a) (Fin c) ℝ) (B : Matrix (Fin a) (Fin d) ℝ)
    (C : Matrix (Fin b) (Fin c) ℝ) (D : Matrix (Fin b) (Fin d) ℝ) :
    r • blk A B C D = blk (r • A) (r • B) (r • C) (r • D) := by
  ext i j
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> simp

lemma castAdd_ne_natAdd' (i : Fin a) (j : Fin b) : Fin.castAdd b i ≠ Fin.natAdd a j := by
  intro h
  have := congrArg Fin.val h
  simp at this
  omega

lemma one_eq_blk : (1 : Matrix (Fin (a + b)) (Fin (a + b)) ℝ) = blk 1 0 0 1 := by
  ext i j
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j
  · simp [Matrix.one_apply]
  · simp [Matrix.one_apply_ne (castAdd_ne_natAdd' i j)]
  · simp [Matrix.one_apply_ne (castAdd_ne_natAdd' j i).symm]
  · simp [Matrix.one_apply]

lemma ones_transpose_mul_ones : (ones a b)ᵀ * ones a c = (a : ℝ) • ones b c := by
  ext i j; simp [ones, Matrix.mul_apply]

/-- `(∑ vᵢ)² ≤ m ∑ vᵢ²`. -/
lemma sum_sq_le_card (v : Fin a → ℝ) : (∑ i, v i) ^ 2 ≤ a * vnorm v ^ 2 := by
  rw [vnorm_sq]
  have := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (fun _ => (1 : ℝ)) v
  simpa [dotProduct, sq] using this

end Cholesky
