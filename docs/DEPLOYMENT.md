# Deployment

Target owner: `caoyhCAS`. Proposed repository: `sheep-reproductive-transcriptomics`.
The source package is version `0.1.0` and includes a GitHub Actions workflow.
It contains no raw sequence reads, credentials, or bundled third-party tools.

If the repository does not yet exist, create an empty repository under that
owner before uploading. Existing repositories for the earlier WGS studies are
different projects and must not be overwritten with this transcriptomics code.

After upload, verify that the default branch contains `README.md`, `config/`,
`scripts/`, `docs/`, `envs/`, `tests/`, and `.github/workflows/ci.yml`. The
`Reconstruction checks` workflow runs the Python tests, base-R input contracts,
synthetic network demo, and the default command planner. The run should finish
successfully; inspect any failure before claiming deployment verified.

CI does not install or execute the historical bioinformatics/statistics stack.
Successful CI verifies helper behavior and contracts, not reproduction of the
paper. Real analyses require the externally obtained reference files, raw reads,
reviewed conflict resolutions and exact legacy dependencies described in the
other documents.

For a manual upload from an extracted source package, one possible route is:

```bash
git init -b main
git add .
git commit -m "Add method reconstruction for Heredity s41437-018-0090-1"
git remote add origin https://github.com/caoyhCAS/sheep-reproductive-transcriptomics.git
git push -u origin main
```

These commands assume the target repository exists and is empty, and that Git
authentication and author identity are configured on your own machine. For a
nonempty target, inspect its history and integrate through a new branch; do not
force-push over existing work.
