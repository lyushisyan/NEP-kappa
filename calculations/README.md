# Material-specific calculations

This is the local research workspace. Everything below it is ignored by Git
except this guide, so raw calculations and personal cluster settings are not
included in a normal source-code commit. Existing files remain on disk.

Recommended layout for each material:

```text
calculations/<material>/
├── README.md          # scientific scope and provenance
├── inputs/            # production YAML/CONTROL/POSCAR inputs (optional)
├── scripts/           # Slurm and post-processing scripts (optional)
├── results/           # generated runs; normally ignored by Git
├── data/              # compact curated tables needed for reproduction
└── figures/           # retained publication figures and metadata
```

Existing historical projects are kept in their current material directories to
avoid breaking provenance. New work should follow the layout above. Reusable
minimal inputs belong in `examples/`; reusable model files belong in
`potentials/<material>/`; package code belongs only in `src/nepkappa/`.

`example-runs/` contains outputs produced by the public examples.
`legacy-examples/` preserves variants removed from the public example catalog
and original site-specific templates. They are local archives, not maintained
public inputs. Small reference summaries required by the automated tests are
copied into `benchmarks/data/`; tests must not depend on this directory.
