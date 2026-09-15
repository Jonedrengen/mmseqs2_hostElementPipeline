# Host Element Pipeline v2: Next Steps

## Recommended Order

1. **Rewrite the Slurm worker for the ten-field manifest**
	- Parse `sample_name`, `trimmed_fasta`, query and reference database prefixes, output paths, coverage modes, and maximum sequence lengths.
	- Keep `SLURM_ARRAY_TASK_ID` supplied by Slurm and use the submitted chunk ID to select exactly one manifest row.

2. **Create one query database per worker**
	- Validate that the trimmed FASTA exists.
	- Create the query database at the manifest-provided prefix before starting searches.
	- Ensure the worker owns the complete lifecycle for one isolate.

3. **Add worker validation and failure handling**
	- Require `-p`, `-e`, `-m`, and `-c` arguments.
	- Validate the manifest is readable and `SLURM_ARRAY_TASK_ID` is set.
	- Fail when no matching manifest row exists.
	- Validate the reference database and required input paths.
	- Stop immediately when query database creation or a search fails.

4. **Make Slurm submission arguments explicit**
	- Pass `pipeline_dir` and `conda_env_prefix` into `start_slurm_runners` rather than relying on globals.
	- Keep the explicit worker script path in the `sbatch` command.

5. **Separate array sizing from concurrency**
	- Replace the hard-coded manifest `max_array_size` with an explicit configuration or function argument.
	- Use the same array-size definition when calculating submission ranges.
	- Keep `max_jobs_per_array` as the array chunk-size limit only.

6. **Use exact chunk lookups**
	- Replace prefix-sensitive `grep` lookups with exact-field `awk` matching for chunk IDs.
	- Confirm that every manifest row belongs to one submitted array range without gaps or duplicates.

7. **Run focused validation before submitting real arrays**
	- Run `bash -n` on the modified scripts.
	- Use a mocked `mmseqs` executable to run one worker against one manifest row.
	- Verify that one query database is created and all configured coverage/length combinations are invoked.
	- Test multiple isolates and a chunk boundary.
	- Submit a small real Slurm array once the local mock test passes.

8. **Defer compiler work**
	- Do not implement `slurm_compiler_worker.sh` until isolate-level runner execution and result layout are verified.
	- Add automatic result combination only after the runner output contract is stable.
