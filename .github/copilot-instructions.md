## Companion Database Repo - Copilot Instructions

## Scope

This repository contains:

- one converter script per external source database
- GitHub workflows to fetch upstream source data from remotes
- generated output artifacts for ResPro import
- validation checks that fail fast on schema or formatting problems

Do not implement profiling logic here. Only fetch, transform, validate, and publish curated
database artifacts.

## Required outputs per external source database

For each external database, produce files in a dedicated folder:

- `rules.tsv` (always required)
- `formula-rules.tsv` (only when combination rules apply)
- `metadata.json` (always required)

Suggested layout:

```text
databases/
	<source_name>/
		scripts/
			convert.py
		output/
			rules.tsv
			formula-rules.tsv   # optional
			metadata.json
```

## ResPro TSV contract

`rules.tsv` must follow the ResPro rules schema used by `respro init --rules`.

Required columns in `rules.tsv`:

- `gene`
- `reference_identifier`
- `position`
- `reference`
- `mutation`
- `antiviral`

Optional columns in `rules.tsv`:

- `group_id`
- `member_id`
- `phenotype`
- `clinical_phenotype`
- `ic50` or `ic_50`
- `fold_ic50` or `fold_ic_50`
- `publication`
- `source`
- `comment`

When `group_id` is used and formula logic is required, emit `formula-rules.tsv` and provide
`member_id` for grouped atomic rows.

`formula-rules.tsv` must follow the ResPro formula schema used by
`respro init --formula-rules`.

Required columns in `formula-rules.tsv`:

- `group_id`
- `antiviral`
- `expression`

Optional columns in `formula-rules.tsv`:

- `label`
- `phenotype`
- `clinical_phenotype`
- `ic50` or `ic_50`
- `fold_ic50` or `fold_ic_50`
- `publication`
- `source`
- `comment`

Expression rules:

- use only `AND`, `OR`, `NOT`, `XOR`
- atomic tokens must match `member_id` from `rules.tsv`
- parentheses are allowed and recommended when precedence is non-trivial

Mutation notation and normalization must stay compatible with ResPro parser expectations.

## Metadata JSON contract

Each generated TSV set must include a `metadata.json` that can be passed to `--metadata`.

Use this fixed top-level object shape:

```json
{
	"maintainers": ["Name or Team"],
	"contact": "email-or-url",
	"publication_pmid": "12345678",
	"website": "https://...",
	"description": "Short description of source and curation scope",
	"maintainer_update": "YYYY-MM-DD",
	"license": "SPDX-or-source-license",
	"tsv_checksum": "sha256:<checksum>"
}
```

Notes:

- `tsv_checksum` must be reproducible and computed from generated TSV artifact content.
- if no PMID exists, use an empty string.
- keep keys stable and do not add ad-hoc fields.

## Converter script rules

For each external source database:

- keep one main converter entrypoint in `databases/<source_name>/scripts/`
- map upstream fields explicitly to ResPro columns; do not infer silently
- enforce deterministic output ordering (stable sort)
- fail fast on unknown mutation syntax or missing required fields
- do not add backward-compatibility layers unless explicitly requested
- avoid network calls inside conversion logic; fetching happens in workflows

If a source does not support combination rules, do not create placeholder
`formula-rules.tsv`.

## GitHub workflow rules

Each external source database must be fetchable from remote via GitHub Actions.

Minimum workflow behavior:

- trigger on schedule and manual dispatch (`workflow_dispatch`)
- fetch upstream source data from remote URL, release asset, or git ref
- run converter script
- generate `rules.tsv`, optional `formula-rules.tsv`, and `metadata.json`
- validate artifact headers and required columns
- compute and inject `tsv_checksum` into `metadata.json`
- upload artifacts (and optionally open PR with updated outputs)

Workflow quality rules:

- pin action versions
- use explicit Python version
- treat warnings about required-column loss as hard failures
- produce clear logs that include upstream source version/date

## Validation checklist

Before merging:

- `rules.tsv` has required columns and no empty required cells
- `formula-rules.tsv` exists only when needed and has valid expressions
- `member_id` and `group_id` consistency is validated
- metadata JSON has only allowed keys and valid date/checksum format
- output is deterministic between repeated runs with same input

## Canonical formatting reference

Primary source of truth for TSV rules semantics and the metadata schema is the `docs/` folder in
the [ResistanceProfiler](https://github.com/the-foxlab/ResistanceProfiler) repository.

Do not rely on undocumented assumptions from external databases.
