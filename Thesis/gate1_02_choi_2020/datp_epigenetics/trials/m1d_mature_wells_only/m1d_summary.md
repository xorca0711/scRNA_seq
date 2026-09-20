# Trial M1d: two rules that could not both hold

**Reading: NOT COMPUTABLE in all eight mature wells**, and for the first
time in this folder the reason was neither biology nor a bad threshold.

M1d removed the neonatal wells from the design, which was right and which
made U7 pass at 3.3-fold enrichment. Then every well refused at the drop
gate, at 27 to 47 per cent, because two of its own frozen rules could not
both be satisfied:

* R5b set the fragment budget to the largest value retaining 60 per cent
  of every well, which guarantees up to 40 per cent of a well falls below it.
* The drop gate inherited from M1b refused any well losing more than 20
  per cent of a group.

A budget defined to drop forty per cent and a gate forbidding more than
twenty cannot both hold. This is the second instance in this repository of
the sub-shape the C12 scope conflict was: two rules written separately,
each defensible, that are not composable.

|  | well | role | n_labelled | reading | why |
|---|---|---|---|---|---|
| 0 | wildtype_PBS | negative control | 18 | not computable | fewer than 50 labelled cells or too small a sham pool |
| 1 | wildtype_SeV | test | 80 | not computable | R5b dropped 47% of a group at the ATAC budget |
| 2 | AP1mut_PBS | negative control | 18 | not computable | fewer than 50 labelled cells or too small a sham pool |
| 3 | AP1mut_SeV | test | 55 | not computable | R5b dropped 34% of a group at the ATAC budget |
| 4 | 7wk_control | negative control | 55 | not computable | R5b dropped 36% of a group at the ATAC budget |
| 5 | 7wk_Cebpa_mutant | negative control | 140 | not computable | R5b dropped 38% of a group at the ATAC budget |
| 6 | SeV_control | test | 44 | not computable | fewer than 50 labelled cells or too small a sham pool |
| 7 | SeV_Cebpa_mutant | test | 400 | not computable | R5b dropped 27% of a group at the ATAC budget |

## Why the fix is better than either rule was

The absolute drop was never the thing worth guarding. If both groups lose
their shallowest fifth, the comparison between them is still fair; what it
loses is generality, and that is a scope limit to report rather than a
defect to refuse on. What corrupts a comparison is DIFFERENTIAL drop. M1e
therefore sets the budget per well from the labelled group and draws the
reference pool only from cells already above it, which makes differential
drop zero by construction rather than merely bounded.
