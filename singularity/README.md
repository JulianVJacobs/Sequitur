**Sequitur Singularity / Apptainer container**

- **Definition file:** `singularity/sequitur.def`
- **Build helper:** `scripts/build_sequitur_sif.sh`
- **SBATCH example:** `singularity/sequitur.sbatch`

Quick notes:

- The definition builds from `alpine:3.18`, installs system build tools via `apk`, Rust, and a few Python packages used by experiments (adjust `pip3 install` as needed).
- The def clones the GitHub repo and builds the Rust release binary `sequitur` and installs it to `/usr/local/bin/sequitur` inside the image.
- To build locally use one of these commands (on a machine with root or apptainer/singularity privileges):

```bash
# with apptainer
apptainer build sequitur.sif singularity/sequitur.def

# or with singularity
singularity build sequitur.sif singularity/sequitur.def
```

- If you want the build to use your current workspace contents rather than cloning GitHub, either:
  - Build a sandbox image and `rsync` your workspace in before converting to a `.sif`, or
  - Edit the `%post` section to use `ADD`/`%files` to copy local files when building with a local builder that supports that.

- Slurm usage: adapt `singularity/sequitur.sbatch` and bind your data dir(s) using `--bind`.

If you want, I can:
- adjust the def to use a different base image (Alpine, Debian), or
- add Python `requirements.txt` support and install exact pinned versions, or
- change the build to produce PyO3 wheels with `maturin` and include Python bindings in the image.
