# `mid_dna_persistence` — chromatin structure as one number

A short note on what this column measures, how it is computed, and what is
and is not known about it.

---

## The idea in one picture

Think of the DNA image inside a nucleus as a landscape, with bright
chromatin as hills.

Now flood it, and slowly lower the water.

- The highest peak appears first, as a single island.
- As the level drops, more islands appear — each is a separate bright
  chromatin domain.
- Keep lowering, and neighbouring islands touch and merge. When two merge,
  we say the *younger* one (the one with the lower peak) has ended.

Every island therefore has a lifetime, measured in brightness: from the level
where it appeared to the level where it merged away. That lifetime is its
**persistence** — not how *bright* the domain is, but how much it *stands out*
from its surroundings.

**`mid_dna_persistence` is the sum of those lifetimes.**

A nucleus with many well-separated chromatin domains accumulates a large sum.
A nucleus whose DNA is gathered into a few large masses accumulates a small
one.

## Why not just count the domains

Counting needs a brightness cutoff: everything above it is a domain,
everything below is not. That makes the answer jumpy. A faint domain sitting
right at the cutoff flips the count by a whole unit, so two nearly identical
nuclei can differ by several counts for no biological reason, and the number
you get depends on a threshold someone chose.

Persistence has no cutoff. A marginal domain contributes a little, a
prominent one contributes a lot, and the total shifts smoothly as the image
changes. It is the continuous version of counting domains.

In the validation described below it also worked better than the count:
separation of condensed from interphase nuclei was 0.75 against 0.63.

## The technical steps

1. **Take the widest plane.** The same mid-section used by the other `mid_`
   columns.
2. **Smooth lightly** (Gaussian, σ ≈ 1 pixel). Without this, every noise
   spike counts as its own island — 500–2800 of them per nucleus, against
   40–220 after smoothing, which is the order of the real domain count.
   Smoothing is a defined operation, not a threshold, so continuity is kept.
3. **Sort every pixel in the nucleus from brightest to darkest**, and add them
   back one at a time. This is the flooding, done exactly rather than in
   steps.
4. **Track which pixels are connected**, using a union-find structure. When a
   newly added pixel touches two existing islands they merge; the one with
   the higher peak survives and the younger one's persistence — its peak
   minus the current level — is recorded.
5. **Normalise and sum.** Each persistence is divided by
   (brightest pixel in the nucleus − field background), making the result
   dimensionless, so it does not change with staining strength, exposure
   time or detector gain. The sum of the normalised values is the reported
   number.

The one island that never merges is left out of the sum. Its persistence is
just the overall contrast of the nucleus, which `mid_dna_cv_corr` already
reports; the islands that *do* merge are the ones carrying information about
internal structure.

No external topology library is used — the implementation is about forty
lines in `nucleus3d/core/quantify.py`.

## How to read the number

| | typical value |
|---|---|
| condensed chromatin (mitotic) | ≈ 1.5 |
| interphase nucleus | ≈ 6.3 |

Lower means fewer, more dominant domains. Higher means many separate
prominent domains.

It is deliberately **independent of chromatin contrast**: measured
correlation with `mid_dna_cv_corr` is r = 0.05. That is the point of having
it. Contrast says how strongly chromatin varies in brightness; persistence
says how that variation is organised into domains. Two nuclei can match on
one and differ on the other, so together they describe more than either
alone.

It also makes two earlier candidate metrics unnecessary: a count of
connected chromatin domains (r = 0.76 with this one) and a skeleton-based
strand-thickness measure (r = 0.64). Neither adds anything once contrast and
persistence are accounted for, so neither was kept.

## What has and has not been shown

**Validation so far.** 268 nuclei from 47 control fields of the cultured-cell
dataset. Each nucleus was scored for cell-cycle phase from its image alone,
with no access to any measured value; 13 came out condensed and 241
interphase, with 14 unscorable. Persistence separates those two groups with
AUC 0.125 — that is, condensed nuclei are reliably *lower* — at
p = 5 × 10⁻⁵ after removing any dependence on chromatin contrast.

**The main limitation is the number of condensed cells.** Thirteen. That is
enough to show the measure works and enough to reject candidates that fail,
but not enough to rank close competitors or to put a confidence interval on
anything.

**The phase labels are model-generated, not ground truth.** They were
produced by scoring images visually. A blind scoring of a subsample by a
person is the planned check and has not been done. Until it is, treat the
numbers above as internally consistent rather than externally validated.

**Most of the condensed nuclei were clipped at the image edge.** The
validation run kept border-touching nuclei (`clear_border=False`, the
package default), and 47% of all nuclei touch an xy edge — but 11 of the 13
condensed ones do, which is more than chance (binomial p = 0.007). Visual
inspection says those are genuine mitotic figures rather than truncation
artifacts: discrete rods with dark gaps between them. Repeating the test
*within* the border-touching nuclei alone, so edge status cannot contribute,
persistence still separates condensed from interphase (AUC 0.229,
p = 0.003; p = 0.003 after removing contrast). The effect is weaker there
than in the full set, so some of the headline separation did come from the
border/interior difference. Restricting to interior nuclei is not possible:
only 2 condensed ones remain. Filtering on `touches_xy_border` for analysis
therefore removes most of the evidence this column was validated against —
worth knowing before relying on it in a filtered table.

**Only control cells were used.** Transcription inhibitors alter chromatin
organisation directly, so in treated samples a cell-cycle signal and a drug
signal cannot be told apart. Whether persistence behaves the same way under
flavopiridol or triptolide is untested.

**Interphase is not sub-staged.** G1, S and G2 cannot be distinguished from
DNA morphology, so this metric helps identify mitotic cells — it does not
order the cell cycle.
