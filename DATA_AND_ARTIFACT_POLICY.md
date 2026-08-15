# Data and Artifact Policy

This repository tracks the reproducible project layer:

- case-builder programs and submission helpers;
- post-processing and plotting scripts;
- research notes, definitions, provenance, and inventory tables;
- curated figures and small reduced tables;
- the two Codex skills used by this project.

Raw HDF5 fields, restart files, archives, videos, and very large long-form
time-series exports are intentionally ignored by Git. They are not source
code and are often many gigabytes. Their source path, case parameters, frame
range, and regeneration method must be recorded in the relevant CSV/README
file before a result is used in a figure.

To check what is excluded locally:

```powershell
git status --ignored --short
```

To reproduce an ignored result, use the script named in its output directory
or rerun the corresponding script in `02_inventory_and_plot_scripts` against
the raw case directory.

