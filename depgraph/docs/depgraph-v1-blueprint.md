# depgraph v1 Blueprint

## Design intent

This document defines the first implementation-ready blueprint for an R package that represents dataset dependency structure as a knowledge graph for leakage-aware biomedical machine learning.

The package is intentionally scoped to dataset dependency graphs, not biological interaction graphs. It is designed to complement `bioLeak`, not duplicate its modeling, preprocessing, resampling, or auditing functions.

Important design choice: the object system uses S3 classes, not S4, per project requirement.

## 1. Package Name Proposals

Final CRAN availability must be checked at submission time. The names below are chosen for clarity, scope fit, and long-term extensibility.

| Candidate | Rationale | Advantages | Disadvantages |
| --- | --- | --- | --- |
| `depgraph` | Short form of dependency graph; matches the core abstraction. | Compact, memorable, broad enough for future extensions. | Generic; could be confused with software dependency graphs. |
| `depGraphR` | Makes the R package identity explicit. | Descriptive, recognizable to R users. | Slightly dated naming style; mixed case is less elegant. |
| `leakgraph` | Emphasizes leakage-relevant graph structure. | Strong alignment with the bioLeak ecosystem. | Too leakage-specific if the package later expands into general provenance analysis. |
| `provgraph` | Focuses on provenance relationships across dataset entities. | Clear scientific framing; extensible beyond ML splitting. | Does not immediately signal dependency constraints. |
| `studygraph` | Highlights study-level provenance and grouping. | Good fit for multi-study biomedical datasets. | Too narrow for subject, batch, and feature provenance. |
| `samplegraph` | Centers the graph on sample-level relationships. | Communicates the unit most users think about first. | Understates higher-level entities like study and feature set. |
| `dataDependR` | Explicitly names data dependencies and has an R-style suffix. | Clear purpose; relatively distinctive. | Less elegant; package name pun may age poorly. |
| `bioDepGraph` | Ties dependency graphs to biomedical datasets. | Strong domain alignment; less ambiguous than `depgraph`. | Slightly longer; future non-biomedical use is less natural. |
| `provsplit` | Connects provenance with cross-validation design. | Strong link to split derivation use case. | Too focused on resampling; weak for general querying and export. |
| `leakdep` | Blends leakage detection with dependency structure. | Short and distinctive. | Name is less intuitive to first-time users. |

### Recommended long-term name

Recommend `depgraph` as the primary package name, with the package title clarifying the domain:

`depgraph: Dataset Dependency Knowledge Graphs for Leakage-Aware Machine Learning`

Rationale:

- It is the cleanest long-term abstraction.
- It leaves room for future extensions beyond leakage audits.
- It matches the central internal object, `dependency_graph`.
- It aligns with the current repository name.

If `depgraph` is unavailable on CRAN, the fallback should be `bioDepGraph`, not `depGraphR`, because it is more distinctive while preserving the scientific scope.

## 2. Mission Statement

`depgraph` represents dataset dependency structure as a typed graph so that biomedical machine learning workflows can encode sample relationships, provenance, and grouping constraints explicitly. The package supports reproducible evaluation design, exposes leakage-relevant dependency patterns, and provides graph-derived constraints that integrate directly with the broader `bioLeak` ecosystem.

## 3. Conceptual Architecture

The package workflow has nine stages.

### 3.1 Metadata ingestion

Input is a sample-level metadata frame plus optional entity tables. At minimum, a dataset must provide `sample_id`. Additional columns such as `subject_id`, `batch_id`, `study_id`, `timepoint`, `assay`, `featureset_id`, and `outcome` are mapped into the canonical schema.

Primary function:

- `ingest_metadata()`

Responsibilities:

- standardize column names
- coerce identifiers to character
- coerce timepoint fields to ordered representations when possible
- record missingness and schema coverage
- attach dataset-level metadata such as source file, cohort name, and version

### 3.2 Node creation

Each supported entity type is converted to a node table with globally unique identifiers. The package generates typed node IDs such as `sample:S001` and `subject:P001` to avoid collisions across entity types.

Primary function:

- `create_nodes()`

Responsibilities:

- deduplicate entities by source identifier
- create canonical node IDs
- attach required and optional attributes
- preserve source identifiers for round-trip export

### 3.3 Edge definition

Edges encode typed provenance and dependency relationships, usually from `Sample` to another entity. Edges are directed in the canonical graph.

Primary function:

- `create_edges()`

Responsibilities:

- generate canonical edge IDs
- validate source and target types
- store relation type and optional metadata
- allow edge-specific annotations such as collection date, provenance confidence, or derivation source

### 3.4 Graph construction

The graph object combines canonical node and edge tables with an internal `igraph` representation. The tables remain the source of truth; the `igraph` object is a computational index.

Primary function:

- `build_dependency_graph()`

Responsibilities:

- merge node and edge sets
- enforce schema consistency
- construct a directed multigraph
- cache summary counts and core indexes

### 3.5 Graph validation

Validation checks structural integrity, semantic consistency, and leakage-relevant cardinality rules.

Primary function:

- `validate_graph()`

Checks:

- unique node IDs
- unique edge IDs
- edge endpoints exist
- allowed source and target node types for each relation
- one-to-one or many-to-one constraints where declared
- missing required attributes by node type
- inconsistent mappings such as one sample linked to multiple subjects unless explicitly allowed
- invalid time ordering metadata

### 3.6 Dependency querying

The graph is queried either directly through neighbor and path operations or through sample-level projections over selected relations.

Primary functions:

- `query_neighbors()`
- `query_paths()`
- `detect_shared_dependencies()`
- `detect_dependency_components()`

Responsibilities:

- find direct provenance neighbors
- find transitive dependency paths
- identify shared subject, batch, study, or feature provenance
- compute connected components among samples under selected dependency rules

### 3.7 Split constraint derivation

The graph is converted into grouping and blocking constraints that `bioLeak` can consume during resampling.

Primary function:

- `derive_split_constraints()`

Responsibilities:

- derive subject-grouped sample partitions
- derive batch-blocked partitions
- derive study leave-out units
- derive temporal ordering rules
- derive composite dependency closures across multiple relation types

### 3.8 Graph visualization

Visualization emphasizes leakage-relevant structure rather than full biological network exploration.

Primary function:

- `visualize_graph()`

Modes:

- full typed graph
- sample-centered ego graph
- projected sample dependency graph
- study or batch subgraph

### 3.9 Graph export

The graph and derived constraints can be exported for audit trails, supplementary materials, or use in external tools.

Primary function:

- `export_graph()`

Supported outputs:

- GraphML
- JSON
- node CSV
- edge CSV
- constraint CSV

## 4. Node Schema

Version 1 supports eight node types. `Sample` is required in every graph. All other node types are optional at the graph level but supported by the schema whenever the corresponding metadata exist.

### Common canonical node fields

Every node row contains:

- `node_id`: globally unique typed identifier, for example `sample:S001`
- `node_type`: one of the supported node classes
- `node_key`: original domain identifier, for example `S001`
- `label`: human-readable label
- `attrs`: named list of type-specific attributes

### 4.1 Sample

Purpose:

- primary observational unit used in model training and evaluation

Required attributes:

- `sample_id`

Optional attributes:

- `subject_id`
- `batch_id`
- `study_id`
- `timepoint_id`
- `assay_id`
- `featureset_id`
- `outcome_id`
- `sample_role`
- `collection_date`
- `source_file`

Example identifiers:

- `S001`
- `TCGA_BRCA_0001`
- `PREDICT.17.T0`

### 4.2 Subject

Purpose:

- entity representing the individual, specimen donor, patient, animal, or donor-derived unit from which one or more samples originate

Required attributes:

- `subject_id`

Optional attributes:

- `subject_type`
- `species`
- `sex`
- `site`
- `enrollment_study`

Example identifiers:

- `P001`
- `PT_0045`
- `Mouse_12`

### 4.3 Batch

Purpose:

- processing, run, plate, center, or sequencing batch that can induce technical dependence

Required attributes:

- `batch_id`

Optional attributes:

- `batch_type`
- `platform`
- `processing_center`
- `run_date`
- `operator`

Example identifiers:

- `B01`
- `SEQ_RUN_2024_07`
- `PLATE_A3`

### 4.4 Study

Purpose:

- study, cohort, trial, project, or collection source that defines provenance boundaries

Required attributes:

- `study_id`

Optional attributes:

- `study_name`
- `study_type`
- `institution`
- `country`
- `study_phase`

Example identifiers:

- `STUDY_A`
- `TCGA`
- `GSE10072`

### 4.5 Timepoint

Purpose:

- ordered collection event or longitudinal measurement point used for time-aware evaluation

Required attributes:

- `timepoint_id`

Optional attributes:

- `time_index`
- `time_value`
- `time_unit`
- `collection_datetime`
- `visit_label`

Example identifiers:

- `T0`
- `Baseline`
- `Day_30`

### 4.6 Assay

Purpose:

- measurement modality or assay platform used to generate the sample's data representation

Required attributes:

- `assay_id`

Optional attributes:

- `modality`
- `platform`
- `vendor`
- `protocol_version`

Example identifiers:

- `RNAseq`
- `Proteomics_LCMS`
- `ClinicalTable`

### 4.7 FeatureSet

Purpose:

- feature definition or processed representation consumed by modeling, including modality-specific matrices, curated panels, or derived feature pipelines

Required attributes:

- `featureset_id`

Optional attributes:

- `featureset_name`
- `modality`
- `feature_count`
- `derivation_scope`
- `derivation_version`
- `generated_from`

Recommended `derivation_scope` values:

- `external`
- `per_training_split`
- `per_dataset`

Example identifiers:

- `FS_RNA_VST_V1`
- `FS_Clinical_20`
- `FS_MultiOmic_Union_V2`

### 4.8 Outcome

Purpose:

- prediction target or observed outcome attached to a sample or subject

Required attributes:

- `outcome_id`

Optional attributes:

- `outcome_name`
- `outcome_type`
- `outcome_value`
- `outcome_scale`
- `observation_level`

Recommended `observation_level` values:

- `sample`
- `subject`

Example identifiers:

- `DX_CASE`
- `RESPONSE_CR`
- `SURVIVAL_1Y`

## 5. Edge Schema

The canonical graph is directed. For dependency detection, the package can project selected relations onto an undirected sample graph.

### Common canonical edge fields

Every edge row contains:

- `edge_id`: globally unique edge identifier
- `from`: source `node_id`
- `to`: target `node_id`
- `edge_type`: typed relation name
- `attrs`: named list of edge-specific attributes

### Core edge types

| Edge type | Source | Target | Meaning |
| --- | --- | --- | --- |
| `sample_belongs_to_subject` | `Sample` | `Subject` | Sample originates from a subject or donor-level unit. |
| `sample_processed_in_batch` | `Sample` | `Batch` | Sample was processed, measured, or normalized within a technical batch. |
| `sample_from_study` | `Sample` | `Study` | Sample belongs to a study, cohort, project, or collection source. |
| `sample_collected_at_timepoint` | `Sample` | `Timepoint` | Sample was collected or observed at a specific longitudinal point. |
| `sample_measured_by_assay` | `Sample` | `Assay` | Sample is represented by a specific assay or modality. |
| `sample_uses_featureset` | `Sample` | `FeatureSet` | Sample is represented using a specific feature definition or processed matrix. |
| `sample_has_outcome` | `Sample` | `Outcome` | Sample has an observed target value. |
| `subject_has_outcome` | `Subject` | `Outcome` | Subject-level target or endpoint attached at the donor level. |

### Optional v1 derived edge types

| Edge type | Source | Target | Meaning |
| --- | --- | --- | --- |
| `timepoint_precedes` | `Timepoint` | `Timepoint` | Ordered relation used for temporal split derivation when explicit time order is available. |
| `featureset_generated_from_study` | `FeatureSet` | `Study` | Feature set was derived using samples from a study or collection source. |
| `featureset_generated_from_batch` | `FeatureSet` | `Batch` | Feature set was derived using a batch-specific processing context. |

The derived edge types are not required for all graphs, but they materially improve leakage-path auditing for feature engineering provenance.

## 6. Graph Object System

Version 1 uses S3 classes. The objects are list-based wrappers over canonical tables and an `igraph` index.

### 6.1 `graph_node_set`

Structure:

```r
list(
  data = data.frame(
    node_id = character(),
    node_type = character(),
    node_key = character(),
    label = character(),
    attrs = I(list())
  ),
  schema_version = "0.1.0",
  source = list()
)
```

Purpose:

- container for validated nodes of one or more supported types

S3 methods:

- `print.graph_node_set()`
- `summary.graph_node_set()`
- `as.data.frame.graph_node_set()`

### 6.2 `graph_edge_set`

Structure:

```r
list(
  data = data.frame(
    edge_id = character(),
    from = character(),
    to = character(),
    edge_type = character(),
    attrs = I(list())
  ),
  schema_version = "0.1.0",
  source = list()
)
```

Purpose:

- container for validated typed relations

S3 methods:

- `print.graph_edge_set()`
- `summary.graph_edge_set()`
- `as.data.frame.graph_edge_set()`

### 6.3 `dependency_graph`

Structure:

```r
list(
  nodes = graph_node_set,
  edges = graph_edge_set,
  graph = igraph::igraph(),
  metadata = list(
    graph_name = NULL,
    dataset_name = NULL,
    created_at = Sys.time(),
    schema_version = "0.1.0"
  ),
  caches = list(
    sample_ids = character(),
    node_index = integer(),
    edge_index = integer()
  )
)
```

Purpose:

- canonical graph object passed through all downstream operations

S3 methods:

- `print.dependency_graph()`
- `summary.dependency_graph()`
- `plot.dependency_graph()`
- `as.igraph.dependency_graph()`

### 6.4 `graph_query_result`

Structure:

```r
list(
  query = character(1),
  params = list(),
  nodes = data.frame(),
  edges = data.frame(),
  table = data.frame(),
  metadata = list()
)
```

Purpose:

- standardized return object for graph traversal and dependency query functions

S3 methods:

- `print.graph_query_result()`
- `summary.graph_query_result()`
- `as.data.frame.graph_query_result()`

### 6.5 `dependency_constraint`

Structure:

```r
list(
  constraint_id = character(1),
  relation_types = character(),
  sample_map = data.frame(
    sample_id = character(),
    group_id = character(),
    reason = character()
  ),
  transitive = TRUE,
  metadata = list()
)
```

Purpose:

- generic grouping or exclusion rule derived from dependency structure

S3 methods:

- `print.dependency_constraint()`
- `summary.dependency_constraint()`
- `as.data.frame.dependency_constraint()`

### 6.6 `split_constraint`

Structure:

```r
list(
  strategy = character(1),
  sample_map = data.frame(
    sample_id = character(),
    split_unit = character(),
    order_index = integer(),
    block_id = character()
  ),
  recommended_bioleak_args = list(),
  metadata = list()
)
```

Purpose:

- resampling-ready constraint object for use by `bioLeak`

S3 methods:

- `print.split_constraint()`
- `summary.split_constraint()`
- `as.data.frame.split_constraint()`

### 6.7 `leakage_constraint`

Structure:

```r
list(
  issue_type = character(1),
  severity = character(1),
  affected_samples = character(),
  evidence = data.frame(),
  recommendation = character(1),
  metadata = list()
)
```

Purpose:

- warning-oriented object describing graph-derived leakage risks

S3 methods:

- `print.leakage_constraint()`
- `summary.leakage_constraint()`
- `as.data.frame.leakage_constraint()`

## 7. Core API

The user-facing API should remain small and explicit. Constructors accept plain `data.frame` inputs and return S3 objects.

### 7.1 `ingest_metadata()`

Purpose:

- standardize user metadata into the canonical column contract

Arguments:

- `data`: sample-level `data.frame`
- `col_map = NULL`: named character vector mapping user columns to canonical names
- `dataset_name = NULL`
- `strict = TRUE`: error on unmapped required columns when `TRUE`

Return value:

- standardized `data.frame` with canonical column names

Example:

```r
meta <- ingest_metadata(raw_meta, col_map = c(sample_id = "sid", subject_id = "pid"))
```

### 7.2 `create_nodes()`

Purpose:

- create a `graph_node_set` from metadata or entity tables

Arguments:

- `data`
- `type`: one of the supported node types
- `id_col`
- `label_col = NULL`
- `attr_cols = NULL`
- `prefix = TRUE`
- `dedupe = TRUE`

Return value:

- `graph_node_set`

Example:

```r
sample_nodes <- create_nodes(
  meta,
  type = "Sample",
  id_col = "sample_id",
  attr_cols = c("subject_id", "batch_id", "study_id", "timepoint_id", "assay_id", "outcome_id")
)
```

### 7.3 `create_edges()`

Purpose:

- create a `graph_edge_set` from relational columns

Arguments:

- `data`
- `from_col`
- `to_col`
- `from_type`
- `to_type`
- `relation`
- `attr_cols = NULL`
- `allow_missing = FALSE`
- `dedupe = TRUE`

Return value:

- `graph_edge_set`

Example:

```r
subject_edges <- create_edges(
  meta,
  from_col = "sample_id",
  to_col = "subject_id",
  from_type = "Sample",
  to_type = "Subject",
  relation = "sample_belongs_to_subject"
)
```

### 7.4 `build_dependency_graph()`

Purpose:

- assemble node and edge sets into a validated `dependency_graph`

Arguments:

- `nodes`: list of `graph_node_set`
- `edges`: list of `graph_edge_set`
- `graph_name = NULL`
- `dataset_name = NULL`
- `validate = TRUE`

Return value:

- `dependency_graph`

Example:

```r
g <- build_dependency_graph(
  nodes = list(sample_nodes, subject_nodes, batch_nodes, study_nodes, time_nodes, assay_nodes, outcome_nodes),
  edges = list(subject_edges, batch_edges, study_edges, time_edges, assay_edges, outcome_edges),
  graph_name = "demo_dependency_graph",
  dataset_name = "DemoStudy"
)
```

### 7.5 `validate_graph()`

Purpose:

- run structural and semantic validation checks

Arguments:

- `graph`
- `checks = c("ids", "references", "cardinality", "schema", "time")`
- `error_on_fail = FALSE`

Return value:

- `graph_validation_result` list with `valid`, `errors`, `warnings`, and `metrics`

Example:

```r
val <- validate_graph(g)
```

### 7.6 `summarize_graph()`

Purpose:

- provide concise graph metrics for auditing and reports

Arguments:

- `graph`
- `by_type = TRUE`

Return value:

- named list or `data.frame` summary

Example:

```r
summarize_graph(g)
```

### 7.7 `query_neighbors()`

Purpose:

- find direct neighbors around one or more nodes

Arguments:

- `graph`
- `node_ids`
- `edge_types = NULL`
- `direction = c("out", "in", "all")`
- `node_types = NULL`

Return value:

- `graph_query_result`

Example:

```r
query_neighbors(g, node_ids = "sample:S001", edge_types = "sample_belongs_to_subject")
```

### 7.8 `query_paths()`

Purpose:

- identify paths between entities to explain dependency chains

Arguments:

- `graph`
- `from`
- `to`
- `edge_types = NULL`
- `mode = c("shortest", "all_shortest")`
- `max_length = NULL`

Return value:

- `graph_query_result`

Example:

```r
query_paths(g, from = "sample:S001", to = "sample:S002")
```

### 7.9 `detect_shared_dependencies()`

Purpose:

- identify samples that share selected dependency types

Arguments:

- `graph`
- `samples = NULL`
- `via = c("Subject", "Batch", "Study", "FeatureSet")`
- `transitive = FALSE`

Return value:

- `graph_query_result`

Example:

```r
detect_shared_dependencies(g, via = c("Subject", "Batch"))
```

### 7.10 `detect_dependency_components()`

Purpose:

- compute connected sample components under selected relation projections

Arguments:

- `graph`
- `via = c("Subject", "Batch", "Study")`
- `transitive = TRUE`
- `min_size = 1`

Return value:

- `dependency_constraint`

Example:

```r
components <- detect_dependency_components(g, via = c("Subject", "Batch"))
```

### 7.11 `derive_split_constraints()`

Purpose:

- convert graph structure into resampling constraints

Arguments:

- `graph`
- `strategy = c("subject_grouped", "batch_blocked", "study_leave_out", "time_ordered", "composite")`
- `via = NULL`
- `transitive = TRUE`
- `time_var = "time_index"`
- `strict = TRUE`

Return value:

- `split_constraint`

Example:

```r
subject_constraint <- derive_split_constraints(g, strategy = "subject_grouped")
batch_constraint <- derive_split_constraints(g, strategy = "batch_blocked")
composite_constraint <- derive_split_constraints(
  g,
  strategy = "composite",
  via = c("Subject", "Batch", "Study")
)
```

### 7.12 `visualize_graph()`

Purpose:

- render the typed graph or a focused projection for inspection

Arguments:

- `graph`
- `focus = c("full", "sample_projection", "ego")`
- `node_types = NULL`
- `edge_types = NULL`
- `layout = "fr"`
- `seed = 1`

Return value:

- invisibly returns plotting data; draws plot as a side effect

Example:

```r
visualize_graph(g, focus = "sample_projection", edge_types = c("sample_belongs_to_subject", "sample_processed_in_batch"))
```

### 7.13 `export_graph()`

Purpose:

- export graph objects and constraints

Arguments:

- `graph`
- `file`
- `format = c("graphml", "json", "nodes_csv", "edges_csv")`

Return value:

- normalized output path

Example:

```r
export_graph(g, "inst/extdata/demo.graphml", format = "graphml")
```

### 7.14 `as_bioleak_split_spec()`

Purpose:

- convert a `split_constraint` into a lightweight object suitable for `bioLeak`

Arguments:

- `constraint`

Return value:

- named list with fields such as `group_ids`, `block_ids`, `order_index`, and `strategy`

Example:

```r
bioleak_spec <- as_bioleak_split_spec(subject_constraint)
```

## 8. Leakage-Relevant Graph Queries

Version 1 should expose focused, domain-specific queries that produce explicit constraints or warnings instead of generic graph traversals only.

### 8.1 Shared-subject detection

Function:

- `detect_shared_dependencies(graph, via = "Subject")`

Use:

- identifies all samples linked to the same subject
- supports grouped cross-validation
- flags repeated-measures leakage risk

Output:

- `dependency_constraint` for grouping

### 8.2 Shared-batch detection

Function:

- `detect_shared_dependencies(graph, via = "Batch")`

Use:

- identifies technical dependence through processing runs, plates, or centers
- supports batch-blocked splits
- surfaces likely confounding between batch and outcome

Output:

- `dependency_constraint`
- optional `leakage_constraint` when batch and outcome align suspiciously

### 8.3 Shared-study detection

Function:

- `detect_shared_dependencies(graph, via = "Study")`

Use:

- identifies provenance overlap across cohorts or recruitment sources
- supports study leave-out evaluation
- flags subjects appearing across multiple studies

Output:

- `dependency_constraint`
- `leakage_constraint` for cross-study subject overlap

### 8.4 Temporal dependency detection

Functions:

- `derive_split_constraints(graph, strategy = "time_ordered")`
- `query_paths()` over `sample_collected_at_timepoint` and `timepoint_precedes`

Use:

- ensures training samples precede test samples
- supports longitudinal evaluation by subject or globally
- detects samples with missing or non-monotone time order

Output:

- `split_constraint` with `order_index`
- `leakage_constraint` when order is incomplete or contradictory

### 8.5 Feature provenance overlap

Function:

- `detect_shared_dependencies(graph, via = "FeatureSet")`

Use:

- identifies samples represented by the same derived feature matrix or globally computed feature definition
- warns when `FeatureSet$derivation_scope == "per_dataset"` and a user attempts fold-based evaluation

Output:

- `leakage_constraint`

### 8.6 Composite dependency closure

Function:

- `detect_dependency_components(graph, via = c("Subject", "Batch", "Study"), transitive = TRUE)`

Use:

- captures indirect leakage pathways, for example sample A shares a subject with B and B shares a batch with C
- provides conservative split units for high-risk studies

Output:

- `dependency_constraint`
- `split_constraint` when passed through `derive_split_constraints(strategy = "composite")`

### 8.7 Outcome-linked provenance warnings

Function:

- internal helper called by `validate_graph()` and `derive_split_constraints()`

Use:

- warns when batches or studies are near-perfect surrogates for outcome labels
- does not replace `bioLeak` auditing, but provides graph-local provenance signals early

Output:

- `leakage_constraint`

## 9. Example End-to-End Workflow

This example shows the intended user experience for a small longitudinal multi-batch dataset.

```r
library(depgraph)

meta <- data.frame(
  sample_id = c("S1", "S2", "S3", "S4", "S5", "S6"),
  subject_id = c("P1", "P1", "P2", "P3", "P3", "P4"),
  batch_id = c("B1", "B2", "B1", "B3", "B3", "B4"),
  study_id = c("ST1", "ST1", "ST1", "ST2", "ST2", "ST2"),
  timepoint_id = c("T0", "T1", "T0", "T0", "T1", "T0"),
  assay_id = c("RNAseq", "RNAseq", "RNAseq", "Proteomics", "Proteomics", "Proteomics"),
  featureset_id = c("FS_RNA_V1", "FS_RNA_V1", "FS_RNA_V1", "FS_PROT_V1", "FS_PROT_V1", "FS_PROT_V1"),
  outcome_id = c("Case", "Case", "Control", "Case", "Case", "Control"),
  stringsAsFactors = FALSE
)

meta <- ingest_metadata(meta, dataset_name = "DemoStudy")

sample_nodes <- create_nodes(
  meta,
  type = "Sample",
  id_col = "sample_id",
  attr_cols = c("subject_id", "batch_id", "study_id", "timepoint_id", "assay_id", "featureset_id", "outcome_id")
)

subject_nodes <- create_nodes(meta, type = "Subject", id_col = "subject_id")
batch_nodes <- create_nodes(meta, type = "Batch", id_col = "batch_id")
study_nodes <- create_nodes(meta, type = "Study", id_col = "study_id")
time_nodes <- create_nodes(
  transform(meta, time_index = c(0L, 1L, 0L, 0L, 1L, 0L)),
  type = "Timepoint",
  id_col = "timepoint_id",
  attr_cols = c("time_index")
)
assay_nodes <- create_nodes(meta, type = "Assay", id_col = "assay_id")
featureset_nodes <- create_nodes(
  transform(meta, derivation_scope = ifelse(grepl("^FS_", featureset_id), "per_dataset", "external")),
  type = "FeatureSet",
  id_col = "featureset_id",
  attr_cols = c("derivation_scope")
)
outcome_nodes <- create_nodes(meta, type = "Outcome", id_col = "outcome_id")

subject_edges <- create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject")
batch_edges <- create_edges(meta, "sample_id", "batch_id", "Sample", "Batch", "sample_processed_in_batch")
study_edges <- create_edges(meta, "sample_id", "study_id", "Sample", "Study", "sample_from_study")
time_edges <- create_edges(meta, "sample_id", "timepoint_id", "Sample", "Timepoint", "sample_collected_at_timepoint")
assay_edges <- create_edges(meta, "sample_id", "assay_id", "Sample", "Assay", "sample_measured_by_assay")
featureset_edges <- create_edges(meta, "sample_id", "featureset_id", "Sample", "FeatureSet", "sample_uses_featureset")
outcome_edges <- create_edges(meta, "sample_id", "outcome_id", "Sample", "Outcome", "sample_has_outcome")

g <- build_dependency_graph(
  nodes = list(
    sample_nodes, subject_nodes, batch_nodes, study_nodes,
    time_nodes, assay_nodes, featureset_nodes, outcome_nodes
  ),
  edges = list(
    subject_edges, batch_edges, study_edges, time_edges,
    assay_edges, featureset_edges, outcome_edges
  ),
  graph_name = "demo_dependency_graph",
  dataset_name = "DemoStudy"
)

validate_graph(g)
summarize_graph(g)

detect_shared_dependencies(g, via = "Subject")
detect_shared_dependencies(g, via = "Batch")

subject_constraint <- derive_split_constraints(g, strategy = "subject_grouped")
study_constraint <- derive_split_constraints(g, strategy = "study_leave_out")
time_constraint <- derive_split_constraints(g, strategy = "time_ordered")
composite_constraint <- derive_split_constraints(
  g,
  strategy = "composite",
  via = c("Subject", "Batch", "Study")
)

visualize_graph(g, focus = "sample_projection")
export_graph(g, file = "demo_dependency_graph.graphml", format = "graphml")

bioleak_spec <- as_bioleak_split_spec(composite_constraint)
```

Expected interpretation:

- `S1` and `S2` must stay grouped because they share subject `P1`.
- `S4` and `S5` must stay grouped because they share subject `P3`.
- `S1` and `S3` share batch `B1`, which may justify a batch-blocked or composite split.
- all `ST1` samples can be held out together for study-level external validation.
- `FS_RNA_V1` and `FS_PROT_V1` are marked as `per_dataset`, so fold-based reuse of those feature sets should produce leakage warnings.

## 10. Implementation Strategy

### 10.1 Internal representation

Use `igraph` internally, but treat canonical node and edge tables as the primary state.

Reasons:

- `igraph` provides efficient traversals, projections, shortest paths, and connected components.
- tables remain easy to validate, print, export, and unit test.
- this avoids the complexity and deployment burden of a graph database.

### 10.2 Recommended package dependencies

Imports:

- `igraph`
- `jsonlite`
- `utils`
- `stats`

Suggests:

- `testthat`
- `knitr`
- `rmarkdown`
- `ggplot2`
- `ggraph`

Keep `dplyr`, `tibble`, and database backends out of Imports for v1 unless implementation pressure clearly justifies them.

### 10.3 Validation logic

Validation should run in three layers.

Layer 1: structural

- node and edge IDs are unique
- endpoints exist
- node and edge types are supported

Layer 2: semantic

- relation signatures are legal
- required attributes exist by node type
- IDs are character and non-empty

Layer 3: leakage-relevant

- samples with multiple subjects
- subjects linked to conflicting subject-level outcomes
- missing or inconsistent time order
- feature sets marked `per_dataset` used in contexts that imply fold reuse
- cross-study subject overlap

### 10.4 Sample projection algorithm

Many leakage constraints are easier to derive on a sample-to-sample projection than on the full typed graph.

Algorithm:

1. subset edges to selected relation types
2. join samples through shared target nodes
3. form an undirected sample graph
4. compute connected components
5. return component IDs as conservative grouping units

This approach supports both direct grouping and transitive closure.

### 10.5 Export formats

- GraphML: use `igraph::write_graph()`
- JSON: export nodes, edges, metadata, and constraints as nested lists via `jsonlite`
- CSV: export flattened node and edge tables plus constraint tables

### 10.6 CRAN-friendly boundaries

Do not require:

- Neo4j
- Spark
- Arrow
- reticulate
- system graph libraries beyond what `igraph` already handles

The package should install and run on standard CRAN platforms with no external services.

### 10.7 Proposed package layout

```text
depgraph/
  DESCRIPTION
  NAMESPACE
  R/
    ingest.R
    constructors.R
    graph-build.R
    validate.R
    query.R
    constraints.R
    visualize.R
    export.R
    bioleak.R
    methods-print.R
    methods-coerce.R
  tests/testthat/
    test-ingest.R
    test-constructors.R
    test-validate.R
    test-query.R
    test-constraints.R
    test-bioleak.R
  vignettes/
    dependency-graphs.Rmd
    leakage-aware-splitting.Rmd
  inst/extdata/
    demo_metadata.csv
```

## 11. Integration with bioLeak

The package should integrate with `bioLeak` through lightweight constraint objects and sample-level vectors, not through duplicated modeling code.

### 11.1 Graph-derived grouping variables

`derive_split_constraints(strategy = "subject_grouped")` should return a `sample_map` that can be converted into a `group_ids` vector for grouped resampling in `bioLeak`.

### 11.2 Batch-aware resampling constraints

`derive_split_constraints(strategy = "batch_blocked")` should return `block_id` values so `bioLeak` can avoid train-test contamination through shared technical runs.

### 11.3 Study-level holdout detection

`derive_split_constraints(strategy = "study_leave_out")` should produce holdout units aligned with external validation and study transferability tests.

### 11.4 Provenance-based leakage warnings

`leakage_constraint` objects should be emitted when the graph indicates obvious evaluation risks, such as:

- subject overlap across studies
- globally derived feature sets reused across folds
- non-monotone longitudinal splits
- batch-outcome entanglement

These warnings are upstream guidance. They do not replace `bioLeak` audits such as permutation-gap analysis or target-leakage scanning.

### 11.5 Graph-informed duplicate detection

The graph can narrow candidate duplicate or near-duplicate comparisons by restricting search to:

- same subject
- same study
- same assay
- same feature set

This makes downstream `bioLeak` duplicate checks more targeted.

### 11.6 Recommended integration contract

`depgraph` should expose:

- `as_bioleak_split_spec(split_constraint)`
- `as.data.frame(split_constraint)`
- `as.data.frame(dependency_constraint)`

`bioLeak` should only need sample-level vectors and metadata tables, not full graph internals.

## 12. Future Extensions

These entities are explicitly out of scope for version 1 but should be supported later without redesigning the core object model:

- `Gene`
- `Protein`
- `Drug`
- `Disease`
- `Pathway`
- `Dataset`

Potential future biological knowledge graph relations:

- `Gene -> regulates -> Protein`
- `Drug -> targets -> Protein`
- `Disease -> associated_with -> Gene`
- `Gene -> participates_in -> Pathway`
- `Gene -> differentially_expressed_in -> Disease`
- `Dataset -> supported_by -> Publication` or evidence resource

Potential future node types for that extension layer:

- `Gene`
- `Protein`
- `Disease`
- `Drug`
- `Pathway`
- `Dataset`

Potential future relation types for that extension layer:

- `associated_with`
- `targets`
- `participates_in`
- `differentially_expressed_in`
- `supported_by`

Potential future upstream data sources:

- `GEO`
- `PubChem`
- `PDB`
- `DisGeNET`
- `KEGG`

Potential future analyses:

- gene-pathway enrichment over graph neighborhoods
- drug-target mapping linked to assay-derived features
- disease-to-dataset provenance overlays
- dataset similarity networks across studies and modalities

The current S3 table-plus-graph design already supports this extension path because new node and edge types can be added without changing the base `dependency_graph` structure.

## 13. Implementation Roadmap

### Phase 1: core graph representation

Deliverables:

- canonical schema definitions
- `ingest_metadata()`
- `create_nodes()`
- `create_edges()`
- `build_dependency_graph()`
- `validate_graph()`
- print and summary methods

Exit criteria:

- can build and validate graphs from sample metadata

### Phase 2: dependency queries

Deliverables:

- `query_neighbors()`
- `query_paths()`
- `detect_shared_dependencies()`
- `detect_dependency_components()`

Exit criteria:

- can explain direct and transitive dependency structure for samples

### Phase 3: split constraint generation

Deliverables:

- `derive_split_constraints()` for subject, batch, study, time, and composite strategies
- `dependency_constraint`, `split_constraint`, and `leakage_constraint` methods

Exit criteria:

- graph-derived sample maps can drive leakage-aware splitting

### Phase 4: bioLeak integration

Deliverables:

- `as_bioleak_split_spec()`
- interoperability tests using representative `bioLeak` workflows

Exit criteria:

- graph outputs feed directly into `bioLeak` resampling interfaces

### Phase 5: visualization and export

Deliverables:

- `visualize_graph()`
- `export_graph()`
- vignette and example data

Exit criteria:

- graphs can be inspected, shared, and archived in publication-ready workflows

## 14. Final v1 Package Blueprint

The proposed version 1 package should be a focused, CRAN-friendly R package for representing dataset dependency structure as a typed provenance graph. It should:

- treat `Sample` as the central unit
- encode subject, batch, study, timepoint, assay, feature set, and outcome dependencies
- store canonical node and edge tables alongside an internal `igraph` index
- use S3 classes for pragmatic scientific package ergonomics
- derive explicit grouping, blocking, and ordering constraints for `bioLeak`
- emit targeted leakage warnings when provenance structure indicates evaluation risk

The package should not fit models, preprocess features, run resampling, or duplicate `bioLeak` auditing. Its role is to make the dependency structure of a dataset explicit, queryable, auditable, and operational for leakage-aware machine learning.

As a first scientific release, this design is narrow enough for a robust CRAN implementation and broad enough to support later method papers on graph-informed leakage prevention, provenance-aware benchmark design, and cross-study reproducibility in biomedical machine learning.
