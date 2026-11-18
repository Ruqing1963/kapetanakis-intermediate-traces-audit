# Prescribed traces in intermediate extensions: audit and corrected bounds

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24345226.svg)](https://doi.org/10.5281/zenodo.24345226)

Paper, verification code and data for

> Ruqing Chen, *Normal and primitive normal elements with prescribed traces in intermediate extensions:
> counterexamples, corrected main term, and improved thresholds*, 2026.
> [doi:10.5281/zenodo.24345226](https://doi.org/10.5281/zenodo.24345226)

It audits

> A. C. Mazumder, G. Kapetanakis, D. K. Basnet, *Normal and primitive normal elements with prescribed traces in
> intermediate extensions of finite fields*, Finite Fields Appl. 110 (2026) 102745,
> [doi:10.1016/j.ffa.2025.102745](https://doi.org/10.1016/j.ffa.2025.102745), arXiv:2510.14840.

## Findings

Notation: `d_1 < ... < d_k` are divisors of `m`, `D = d_1 + ... + d_k`, `g = lcm(x^{d_i} - 1)`, `lambda = deg g`.

1. The count of normal elements with prescribed traces (Theorem 3.6) is correct, also without the paper's
   standing hypotheses.
2. The main term of the character sum in Section 4 loses a factor `q^D`. The sufficient condition of
   Theorem 4.1 is printed with exponent `m/2 - lambda - D`; the correct exponent is `m/2 - lambda`.
3. Theorem 5.1 omits the hypothesis that the prescribed traces be normal. As printed it fails for exactly
   `q^lambda - Phi(g) > 0` admissible trace tuples, for every choice of parameters.
4. Two threshold proofs in Theorem 5.1 (`k = 4` and `k = 3`) rest on inequalities that fail for large `m`.
5. Theorem 4.2 (the case `lcm(d_i) < m`) is not proved by the indicated argument. It is proved here when
   `m / lcm(d_i)` is a power of the characteristic and whenever the corrected condition applies; otherwise it
   is open.
6. The reduction of non-coprime degrees to coprime ones by a change of base field is invalid.

Corrected thresholds for pairwise coprime degrees and normal admissible traces:

| k       | printed                              | proved here                         |
|---------|--------------------------------------|-------------------------------------|
| 2       | `q` large, `d_1 >= 6` with exclusions | `q` large, `d_1 >= 3`, `d != (3,4)` |
| 3       | `q >= 2.27e24072855`, `m >= 60`      | `q >= 8003`, all `m`                |
| 4       | `q >= 1334`                          | `q >= 4`                            |
| 5, 6, 7 | `q >= 9, 7, 5`                       | `q >= 4`                            |

No counterexample to an existence statement was found once the prescribed traces are assumed normal.
The cases `k >= 4` with `q in {2, 3}` and `k = 3` with `q < 8003` remain open.

## Contents

| Path | Content |
|------|---------|
| `mazumder_kapetanakis_intermediate_traces_audit.tex`, `.pdf` | the paper |
| `fig_intermediate_traces_audit.pdf` | the figure |
| `Code/verify_intermediate_traces.py` | exhaustive enumeration in `F_{q^m}` and threshold computations |
| `Code/verify_base_change_example.py` | explicit example in `F_{5^12}` (Proposition 5.2(b)) |
| `Code/make_figure.py` | builds the figure from the data |
| `Data/exhaustive.json`, `Data/exhaustive_log.txt` | results of the 36 exhaustive cases (`k = 2`, `q <= 13`, `q^m <= 2^24`) |
| `Data/numeric.json`, `Data/numeric_log.txt` | constants, checks of the printed inequalities, corrected thresholds |
| `Data/base_change_example.json` | the element of `F_{5^12}` and its verified properties |

## Reproducing

Requires Python 3 with `numpy`, `sympy` and `matplotlib` (run with Python 3.12, numpy 1.26, sympy 1.13,
matplotlib 3.8). Run from the `Code` directory; outputs are written to `Data`.

```
python verify_intermediate_traces.py exhaustive small   # fields up to 2^22 elements, under a minute
python verify_intermediate_traces.py exhaustive all     # adds the fields with up to 2^24 elements; a few minutes, several GB of RAM
python verify_intermediate_traces.py numeric            # constants and thresholds, about a minute
python verify_intermediate_traces.py numeric c31        # also recomputes C_31 by sieving to 2^31 (slow)
python verify_base_change_example.py
python make_figure.py
```

Without the `c31` option the stored value `log10 C_31 = 1553069.3804745` is used; it agrees with the value
printed in the audited paper.

## Citation

```bibtex
@misc{chen2026intermediatetraces,
  author = {Chen, Ruqing},
  title  = {Normal and primitive normal elements with prescribed traces in intermediate extensions:
            counterexamples, corrected main term, and improved thresholds},
  year   = {2026},
  doi    = {10.5281/zenodo.24345226},
  url    = {https://doi.org/10.5281/zenodo.24345226}
}
```

## License

The paper, figure and data are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); the code in `Code/` is licensed under the MIT
License. See `LICENSE`.
