"""
Joining the toolbox's two tables: nuclei, and what is inside them.

This module holds the join and nothing else. It lives outside both packages
on purpose. `antenna3d` was written to be extractable and usable on its own,
with masks from anywhere, and `nucleus3d` does not need to know that anything
downstream exists; an import between them would end both properties. The join
is pandas and a handful of consistency checks, so it costs nothing to keep it
here.

    from combine import join_nuclei_and_antennas
    both = join_nuclei_and_antennas("results/nuclei", "results/antennas")
"""

import os

import pandas as pd

#: Columns both tables carry with the same meaning and the same value. Verified
#: row by row, then kept once instead of being suffixed into two.
_SHARED = ("position", "label")


def _read(path, what):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"no {what} at {path}. Run that half of the pipeline first -- "
            f"scripts/run_segmentation.py writes the nucleus table, "
            f"scripts/run_antennas.py the antenna table, and "
            f"scripts/run_combined.py runs both in order.")
    return pd.read_csv(path)


def _check_no_failures(outdir, log):
    """A run that lost fields produces a table that looks complete.

    Both pipelines record a failed field and carry on, which is the right
    behaviour -- one unreadable file must not cost the other twelve. It does
    mean the tables are silent about what is missing, and `failures.csv` is
    the only place it shows.
    """
    fp = os.path.join(outdir, "failures.csv")
    if os.path.exists(fp):
        n = len(pd.read_csv(fp))
        log(f"  WARNING: {fp} exists -- {n} field(s) failed and are absent "
            f"from this join.")


def join_nuclei_and_antennas(nucleus_dir, antenna_dir, out_csv=None, log=print):
    """
    One row per nucleus that both pipelines measured.

    `nucleus_dir` is a nucleus3d output tree, `antenna_dir` an antenna3d one,
    run over the same data. The join key is `nucleus_uid`, which both build
    forward from (file stem, position, label).

    **Fewer rows than either input is correct, not a bug.** The two pipelines
    apply different gates on purpose: nucleus3d keeps small objects
    (`min_volume_um3` defaults to 15, because a slab volume is not a nuclear
    volume) and keeps nuclei touching the image edge, while antenna3d gates at
    100-4000 um3 and drops xy-border-touching nuclei by default, because a
    nucleus cut laterally has no usable denominator. The reconciliation is
    logged rather than hidden.

    An EMPTY overlap raises. A silent zero-row merge is the failure this join
    exists to prevent: it is what a mismatched `nucleus_uid` produces, and it
    looks exactly like a dataset with no antennas in it.

    Returns
    -------
    DataFrame, with the counts also available as `df.attrs["join"]`.
    """
    nuc = _read(os.path.join(nucleus_dir, "nuclei_measurements.csv"),
                "nucleus table")
    ant = _read(os.path.join(antenna_dir, "antenna_nuclei.csv"),
                "antenna table")

    _check_no_failures(nucleus_dir, log)
    _check_no_failures(antenna_dir, log)

    for name, df in (("nucleus", nuc), ("antenna", ant)):
        if "nucleus_uid" not in df.columns:
            raise KeyError(f"the {name} table has no nucleus_uid column")
        dup = df.nucleus_uid.duplicated().sum()
        if dup:
            raise ValueError(
                f"the {name} table has {dup} duplicated nucleus_uid value(s). "
                f"The id is unique per (file stem, position, label); duplicates "
                f"mean two source files share a stem, and the join would "
                f"multiply rows rather than match them.")

    # Suffixed, not collapsed. `file` is a basename with its extension on one
    # side and a stem on the other, `source_path` depends on how each run was
    # invoked, and `condition` spells "no subfolder" differently. Only the
    # columns checked below are safe to treat as one.
    both = nuc.merge(ant, on="nucleus_uid", how="inner",
                     suffixes=("_nucleus", "_antenna"))

    if not len(both):
        raise ValueError(
            f"no nucleus_uid is present in both tables: {len(nuc)} nucleus rows "
            f"and {len(ant)} antenna rows share none.\n"
            f"  a nucleus id: {nuc.nucleus_uid.iloc[0] if len(nuc) else '(none)'}\n"
            f"  an antenna id: {ant.nucleus_uid.iloc[0] if len(ant) else '(none)'}\n"
            f"Either the two runs were over different data, or the antenna run "
            f"segmented its own nuclei instead of reading nucleus3d's labels "
            f"(check that LABELS_DIR pointed at {nucleus_dir}/labels).")

    # The columns that must agree, collapsed back to one after checking they do.
    for col in _SHARED:
        a, b = f"{col}_nucleus", f"{col}_antenna"
        if a in both.columns and b in both.columns:
            bad = both[both[a] != both[b]]
            if len(bad):
                raise ValueError(
                    f"{len(bad)} row(s) disagree on `{col}` between the two "
                    f"tables while sharing a nucleus_uid. The id is not "
                    f"identifying what it claims to; do not use this join. "
                    f"First: {bad.nucleus_uid.iloc[0]}")
            both = both.drop(columns=[b]).rename(columns={a: col})

    front = [c for c in ("nucleus_uid", "position", "label",
                         "file_nucleus", "condition_nucleus") if c in both.columns]
    both = both[front + [c for c in both.columns if c not in front]]

    stats = {"n_nuclei": len(nuc), "n_antenna": len(ant), "n_joined": len(both),
             "n_nucleus_only": len(nuc) - len(both),
             "n_antenna_only": len(ant) - len(both)}
    both.attrs["join"] = stats

    log(f"  joined {stats['n_joined']} nuclei "
        f"({stats['n_nuclei']} in the nucleus table, {stats['n_antenna']} in the "
        f"antenna table)")
    if stats["n_nucleus_only"]:
        log(f"  {stats['n_nucleus_only']} nucleus row(s) have no antenna row. "
            f"Expected: antenna3d gates on volume and drops nuclei touching an "
            f"xy edge, and nucleus3d does neither by default.")
    if stats["n_antenna_only"]:
        log(f"  {stats['n_antenna_only']} antenna row(s) have no nucleus row. "
            f"NOT expected -- the antenna run may have segmented its own nuclei "
            f"rather than reading nucleus3d's labels.")

    if out_csv:
        os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
        both.to_csv(out_csv, index=False)
        log(f"  wrote {out_csv}")

    return both
