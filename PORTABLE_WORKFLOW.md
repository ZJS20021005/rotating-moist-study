# Portable Git Workflow

## Repository layout

The main repository is the project record and analysis layer:

```text
rotating_case_inventory/
  00_latest_program/              DNS source/build artefacts
  01_case_builder_remote_upload/  case creation and submission tools
  02_inventory_and_plot_scripts/  deterministic analysis and plotting
  03_inventory_tables/            case bookkeeping
  04_outputs_and_figures/         curated figures and reduced results
  05_source_reference/            papers, notes, and Codex skill mirrors
```

Large simulation fields remain on the data drive or cluster. Git records how
to find and regenerate them rather than duplicating them into the repository.

## Use on another computer

1. Clone this repository.
2. Open the cloned directory in VS Code.
3. Run:

```powershell
.\install_codex_skills.ps1
```

4. Edit the machine-specific paths in
`01_case_builder_remote_upload/case_config.json` and in any plotting script
that still points to a local data drive.
5. Run the case-builder or plotting scripts from the cloned repository.

The physical parameters and analysis definitions stay in Git. Only machine
paths, SSH configuration, cluster modules, and raw data locations change.

## Update and share changes

```powershell
git pull --rebase
git add .
git commit -m "Update case records and diagnostics"
git push
```

Always inspect `git status` before committing. Never commit SSH keys,
passwords, raw restart files, or cluster credentials.

