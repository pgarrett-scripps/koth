# Feature finding

```bash
# mzML input, results written to ./out/<stem>/
koth_ff data.mzML --output ./out

# Bruker .d directory
koth_ff data.d --output ./out

# Thermo .raw (native reader, default build)
koth_ff data.raw --output ./out

# Gzip-compressed mzML is streamed
koth_ff data.mzML.gz --output ./out

# With a custom config
koth_ff data.mzML --config my_config.toml --output ./out

# Skip the scoring stage (faster)
koth_ff data.mzML --output ./out --no-scoring
```

The output directory `<output>/<input_stem>/` contains the following with the default TSV output:

```
hills.tsv
features.tsv
report.json     # summary statistics (hill and feature counts, score distribution) and the 1/K0 scale used
config.toml     # the config that was used, for reproducibility
```

Set `[output].format = "parquet"` in your TOML to write `hills.parquet` and
`features.parquet` instead. See [output formats](OUTPUTS.md) and the
[configuration reference](CONFIGURATION.md) for the optional MS2 outputs,
preprocessing, and instrument settings.

For Bruker `.d` input, every reported ion mobility is Bruker's
acquisition-calibrated 1/K0 (`[file] bruker_mobility_scale = "calibrated"`,
the default since 0.10.0). `"linear"` restores the timsrust converter used up
to 0.9.0. `koth_align` refuses a batch whose runs used different scales.

Start with [example_config.toml](../example_config.toml). Every setting is
optional; unknown fields are rejected. The resolved configuration is saved
with the results, including CLI overrides.

For multi-run analysis, see [alignment and LFQ](LFQ.md).
