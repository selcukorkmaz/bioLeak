# depgraph

`depgraph` is an R package for representing dataset dependency structure as a typed knowledge graph for leakage-aware biomedical machine learning.

The package is designed for the `bioLeak` ecosystem. It does not fit models, preprocess features, or run resampling itself. Its job is to make dataset structure explicit so that subject overlap, batch effects, study provenance, time ordering, assay provenance, feature provenance, and outcome linkage can be encoded and audited before model evaluation.

## Package status

Current status: early scaffold / alpha.

Implemented now:

- package metadata and installable package structure
- S3 object system for graph containers and constraint objects
- centralized node and edge schema definitions
- metadata ingestion
- typed node creation
- typed edge creation
- dependency graph assembly with `igraph`
- structured structural, semantic, and leakage-relevant validation
- dependency query and traversal functions
- projected sample dependency detection helpers
- split-constraint derivation for subject, batch, study, time, and composite modes
- bioLeak-facing split-spec translation and preflight validation
- structural leakage risk summaries for graph and constraint handoff
- print, summary, table coercion, and plotting methods
- compatibility constructor aliases for object-system integration
- manual pages, demo metadata, and unit tests

Planned next from the v1 blueprint:

- leakage-path detection helpers
- graph export helpers

## Why this package exists

Leakage in biomedical machine learning often comes from the structure of the dataset, not only from mistakes in code. Samples can share:

- the same subject
- the same batch
- the same study
- the same collection timepoint
- the same assay
- the same derived feature set
- the same outcome definition

If those dependencies are not encoded explicitly, train-test separation can look valid while still being scientifically invalid.

`depgraph` exists to make those relationships first-class objects.

## Scope

Version 1 is intentionally limited to dataset dependency graphs. It is not a biological interaction graph package.

Included entity types:

- `Sample`
- `Subject`
- `Batch`
- `Study`
- `Timepoint`
- `Assay`
- `FeatureSet`
- `Outcome`

Not included in v1:

- gene networks
- protein interaction graphs
- pathway graphs
- drug-target graphs

Those are future extensions, not part of the current implementation.

## Core features

### 1. Metadata ingestion

Function:

- `ingest_metadata()`

What it does:

- standardizes sample metadata into a canonical schema
- supports optional column remapping through `col_map`
- coerces ID-like fields to character
- stores `dataset_name` as metadata on the returned table

Why it matters:

- downstream graph construction depends on consistent identifiers and stable column names

### 2. Typed node creation

Function:

- `create_nodes()`

What it does:

- converts ordinary metadata tables into canonical node tables
- creates globally unique typed IDs such as `sample:S1` and `subject:P1`
- stores extra attributes in an `attrs` list-column
- deduplicates repeated entities

Why it matters:

- prevents ID collisions across entity types
- turns flat metadata into graph-ready objects

### 3. Typed edge creation

Function:

- `create_edges()`

What it does:

- builds canonical relation tables from source and target columns
- validates source and target node types for known edge types
- creates stable typed edges such as `sample_belongs_to_subject`
- supports optional edge attributes

Why it matters:

- graph meaning lives in the edge types
- edge typing is what later enables dependency queries and split constraints

### 4. Dependency graph assembly

Functions:

- `build_dependency_graph()`
- `build_depgraph()`
- `as_igraph()`

What it does:

- combines node and edge sets into a single `dependency_graph`
- stores canonical node and edge tables as source of truth
- builds an internal `igraph` index for traversal and graph algorithms
- stores graph metadata and lookup caches

Why it matters:

- gives you one object that can be printed, summarized, validated, plotted, or converted to `igraph`

### 5. Graph validation

Functions:

- `validate_graph()`
- `validate_depgraph()`

What it checks:

- duplicate node IDs
- duplicate edge IDs
- duplicate illegal `(from, to, edge_type)` relations
- missing edge endpoints
- unsupported node types
- illegal source-target signatures for known relations
- invalid or unknown edge types
- single-target relation violations such as sample-to-subject
- invalid `timepoint_precedes` self-loops
- malformed `time_index` values when present
- invalid `FeatureSet$derivation_scope`
- invalid `Outcome$observation_level`
- repeated-subject sample structures
- cross-study subject overlap
- per-dataset feature derivation scope
- graph-local shared feature provenance risks

Why it matters:

- catches structural errors, semantic inconsistencies, and graph-local leakage risk signals before downstream evaluation design

Validation output:

- returns a `depgraph_validation_report`
- classifies issues as `error`, `warning`, or `advisory`
- preserves backward-compatible fields such as `$valid`, `$errors`, `$warnings`, and `$metrics`
- supports explicit validation overrides through `build_dependency_graph(..., validation_overrides = ...)`

### 6. S3 object system

Constructors:

- `graph_node_set()`
- `graph_edge_set()`
- `dependency_graph()`
- `new_depgraph_nodes()`
- `new_depgraph_edges()`
- `new_depgraph()`
- `graph_query_result()`
- `dependency_constraint()`
- `split_constraint()`
- `leakage_constraint()`

What it does:

- provides lightweight, CRAN-friendly, scientific-package style containers
- keeps the package S3-based instead of S4
- gives each object predictable fields and methods

Why it matters:

- keeps the package simple to inspect, extend, and integrate with `bioLeak`

### 7. Inspection and coercion methods

Available methods:

- `print()`
- `summary()`
- `as.data.frame()`
- `plot()` for `dependency_graph`

What they do:

- print concise object summaries
- report node and edge counts by type
- expose canonical tables for export or manual inspection
- plot the current `igraph` representation

Why it matters:

- the package is immediately usable for inspection and graph-level exploration

### 8. Dependency query layer

Functions:

- `query_node_type()`
- `query_edge_type()`
- `query_neighbors()`
- `query_paths()`
- `query_shortest_paths()`
- `detect_dependency_components()`
- `detect_shared_dependencies()`

What it does:

- retrieves typed node and edge subsets from the canonical tables
- traverses graph neighborhoods and dependency paths through the `igraph` backend
- returns standardized `graph_query_result` objects
- detects projected sample dependency components through shared subjects, batches, studies, timepoints, assays, feature sets, or outcomes
- identifies direct sample pairs that share dependency sources such as the same subject, batch, study, or timepoint

Why it matters:

- makes dataset dependency structure directly queryable before any split-constraint logic is derived

### 9. Documentation, test, and demo data

Included assets:

- manual pages in `man/`
- smoke test in `tests/testthat/`
- demo metadata in `inst/extdata/demo_metadata.csv`
- design blueprint in `docs/depgraph-v1-blueprint.md`

Why it matters:

- the scaffold is already usable as a starting point for real package development

### 10. Split-constraint generation

Functions:

- `derive_split_constraints()`
- `grouping_vector()`

Supported modes:

- `subject`
- `batch`
- `study`
- `time`
- `composite`

What it does:

- derives sample-level grouping assignments from dependency structure
- creates reusable `split_constraint` objects
- supports strict composite closure over projected sample dependency graphs
- supports deterministic rule-based composite grouping with dependency priorities
- exposes a named grouping vector for downstream cross-validation interfaces

Why it matters:

- turns dependency graphs into operational grouping specifications that later `bioLeak` workflows can consume

How constraints are derived:

- `subject`: samples sharing the same `Subject` node through `sample_belongs_to_subject` receive the same `group_id`
- `batch`: samples sharing the same `Batch` node through `sample_processed_in_batch` receive the same `group_id`; missing batch links are retained as singleton unlinked groups with warnings
- `study`: samples sharing the same `Study` node through `sample_from_study` receive the same `group_id`
- `time`: samples are grouped by `Timepoint`; ordering is derived from `Timepoint$time_index` when present, otherwise from `timepoint_precedes` edges when the timepoint graph is a DAG
- `composite` + `strict`: selected dependency relations are projected onto a sample graph and connected components become the split groups
- `composite` + `rule_based`: each sample is assigned by the highest-priority available dependency source; lower-priority available sources are preserved in the explanation field

### 11. bioLeak integration layer

Functions:

- `as_bioleak_split_spec()`
- `validate_bioleak_split_spec()`
- `summarize_leakage_risks()`

What it does:

- translates `split_constraint` objects into a stable `bioleak_split_spec`
- preserves grouping, blocking, and time-ordering columns in a canonical `sample_data` table
- performs preflight structural checks before any downstream `bioLeak` resampling is attempted
- summarizes graph validation issues, constraint warnings, and split-spec readiness in a human-readable form

Why it matters:

- gives `bioLeak` a clean handoff object without making `depgraph` responsible for model fitting or fold generation

## Supported node types

Each node has the canonical fields:

- `node_id`
- `node_type`
- `node_key`
- `label`
- `attrs`

Supported node types and intended meaning:

| Node type | Meaning |
| --- | --- |
| `Sample` | Primary observational unit used in model evaluation. |
| `Subject` | Donor, patient, animal, or donor-derived unit linked to one or more samples. |
| `Batch` | Processing or measurement batch that can induce technical dependence. |
| `Study` | Cohort, trial, or dataset provenance boundary. |
| `Timepoint` | Ordered longitudinal collection event. |
| `Assay` | Measurement modality or assay platform. |
| `FeatureSet` | Derived feature representation used for modeling. |
| `Outcome` | Prediction target or observed endpoint. |

## Supported edge types

Each edge has the canonical fields:

- `edge_id`
- `from`
- `to`
- `edge_type`
- `attrs`

Core supported relations:

| Edge type | Meaning |
| --- | --- |
| `sample_belongs_to_subject` | Sample originates from a subject. |
| `sample_processed_in_batch` | Sample was processed in a batch. |
| `sample_from_study` | Sample belongs to a study or cohort. |
| `sample_collected_at_timepoint` | Sample was observed at a timepoint. |
| `sample_measured_by_assay` | Sample was generated by an assay. |
| `sample_uses_featureset` | Sample is represented by a feature set. |
| `sample_has_outcome` | Sample has an outcome value. |
| `subject_has_outcome` | Subject has a subject-level outcome. |
| `timepoint_precedes` | One timepoint precedes another. |
| `featureset_generated_from_study` | Feature set was derived in a study context. |
| `featureset_generated_from_batch` | Feature set was derived in a batch context. |

## Object model

The package currently centers on seven S3 classes.

| Class | Purpose |
| --- | --- |
| `graph_node_set` | Canonical node table wrapper. |
| `graph_edge_set` | Canonical edge table wrapper. |
| `dependency_graph` | Main graph object containing tables plus an `igraph` index. |
| `graph_query_result` | Standardized return type for query and dependency-detection APIs. |
| `dependency_constraint` | Generic dependency grouping object. |
| `split_constraint` | Resampling-ready constraint object for future split derivation. |
| `leakage_constraint` | Warning-oriented object for future leakage detection. |

## Installation

From the package root in R:

```r
install.packages(".", repos = NULL, type = "source")
```

Then load it with:

```r
library(depgraph)
```

## Quick start

```r
library(depgraph)

meta <- data.frame(
  sample_id = c("S1", "S2", "S3", "S4"),
  subject_id = c("P1", "P1", "P2", "P3"),
  batch_id = c("B1", "B2", "B1", "B3"),
  study_id = c("ST1", "ST1", "ST1", "ST2"),
  stringsAsFactors = FALSE
)

meta <- ingest_metadata(meta, dataset_name = "DemoStudy")

sample_nodes <- create_nodes(
  meta,
  type = "Sample",
  id_col = "sample_id",
  attr_cols = c("subject_id", "batch_id", "study_id")
)

subject_nodes <- create_nodes(meta, type = "Subject", id_col = "subject_id")
batch_nodes <- create_nodes(meta, type = "Batch", id_col = "batch_id")
study_nodes <- create_nodes(meta, type = "Study", id_col = "study_id")

subject_edges <- create_edges(
  meta,
  from_col = "sample_id",
  to_col = "subject_id",
  from_type = "Sample",
  to_type = "Subject",
  relation = "sample_belongs_to_subject"
)

batch_edges <- create_edges(
  meta,
  from_col = "sample_id",
  to_col = "batch_id",
  from_type = "Sample",
  to_type = "Batch",
  relation = "sample_processed_in_batch"
)

study_edges <- create_edges(
  meta,
  from_col = "sample_id",
  to_col = "study_id",
  from_type = "Sample",
  to_type = "Study",
  relation = "sample_from_study"
)

g <- build_dependency_graph(
  nodes = list(sample_nodes, subject_nodes, batch_nodes, study_nodes),
  edges = list(subject_edges, batch_edges, study_edges),
  graph_name = "demo_dependency_graph",
  dataset_name = "DemoStudy"
)

print(g)
summary(g)
validate_graph(g)

node_table <- as.data.frame(g$nodes)
edge_table <- as.data.frame(g$edges)
ig <- as_igraph(g)

plot(g)
```

## Example output interpretation

In the quick-start example:

- `S1` and `S2` share the same subject `P1`
- `S1` and `S3` share the same batch `B1`
- `S1`, `S2`, and `S3` come from the same study `ST1`

Those are precisely the dependency patterns that later split-constraint logic will need to account for.

## What is implemented vs planned

Current implementation:

- canonical schema constants
- node and edge constructors
- dependency graph constructor
- multi-level validation framework with structured reports
- dependency query and path traversal APIs
- projected sample dependency detection
- split-constraint derivation and grouping-vector extraction
- bioLeak split-spec translation, preflight validation, and leakage summaries
- `igraph` conversion
- S3 methods
- package docs and unit tests

Not implemented yet:

- `visualize_graph()` wrapper beyond the base `plot()` method
- `export_graph()`

The core query, constraint, and `bioLeak` translation layers now exist, and the next package milestones are export helpers and richer end-user examples.

## Relationship to bioLeak

`depgraph` is intended to be the structural layer under leakage-aware modeling workflows.

`depgraph` should answer:

- which samples are connected by subject provenance
- which samples share a technical batch
- which samples belong to the same study
- which feature representations were derived globally
- which dependencies imply grouped or blocked resampling

`bioLeak` should answer:

- how to run guarded preprocessing
- how to run leakage-aware resampling
- how to train and evaluate models
- how to audit leakage inflation

That separation keeps both packages focused.

## Repository layout

| Path | Purpose |
| --- | --- |
| `R/` | Core package implementation. |
| `man/` | Manual pages. |
| `tests/testthat/` | Smoke tests. |
| `inst/extdata/` | Demo metadata. |
| `docs/depgraph-v1-blueprint.md` | Full v1 design and implementation blueprint. |

## Development roadmap

The next milestones are:

- implement graph export helpers
- add vignettes and richer examples

For the full design specification, see [`docs/depgraph-v1-blueprint.md`](docs/depgraph-v1-blueprint.md).
