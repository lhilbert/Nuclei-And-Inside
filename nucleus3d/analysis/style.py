"""
Shared figure constants for the analysis half.

Note what this module does NOT do: it does not call `matplotlib.use()`.
A library that forces a backend takes the choice away from the caller, and
breaks interactive use in a notebook. Matplotlib already falls back to Agg
when no display is available, so headless scripts work without help; the
entry points under `scripts/` set the backend explicitly for themselves.
"""

PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7",
           "#56B4E9", "#D55E00", "#F0E442", "#999999"]
"""
Okabe-Ito, which stays distinguishable under the common forms of
colour-vision deficiency. A red/green pairing does not.
"""

CONTEXT_GREY = "#D9D9D9"
"""For the all-nuclei background layer that gives each panel its context."""

HIGHLIGHT = "#D55E00"
"""Marking picked nuclei in the gallery."""


def sizes(base=9):
    """Font sizes derived from one number, so a figure scales coherently."""
    return dict(base=base, annot=base - 1, tick=base - 2)
