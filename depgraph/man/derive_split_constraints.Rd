\name{derive_split_constraints}
\alias{derive_split_constraints}
\alias{grouping_vector}
\title{Derive Split Constraints from Dependency Graphs}
\usage{
derive_split_constraints(graph,
  mode = c("subject", "batch", "study", "time", "composite"),
  samples = NULL, strategy = c("strict", "rule_based"),
  via = NULL, priority = NULL, include_warnings = TRUE)

grouping_vector(x)
}
\arguments{
  \item{graph}{A \code{dependency_graph}.}
  \item{mode}{Constraint derivation mode.}
  \item{samples}{Optional sample identifiers or sample node IDs used to restrict the returned \code{sample_map}.}
  \item{strategy}{Composite grouping strategy. Ignored for non-composite modes.}
  \item{via}{Optional dependency sources used for composite grouping. May be given as lower-case modes such as \code{"subject"} or node types such as \code{"Subject"}.}
  \item{priority}{Optional priority order used for \code{strategy = "rule_based"}.}
  \item{include_warnings}{Whether to retain human-readable warnings in the returned metadata.}
  \item{x}{A \code{split_constraint}.}
}
\value{
\code{derive_split_constraints()} returns a \code{split_constraint} whose
\code{sample_map} contains grouping assignments and, for time-aware
constraints, ordering metadata. \code{grouping_vector()} returns a named
character vector of \code{group_id} values keyed by \code{sample_id}.
}
\description{
Convert dataset dependency structure into deterministic sample-level grouping
constraints suitable for leakage-aware evaluation design.
}
\details{
Constraint derivation rules:

\itemize{
  \item{\code{mode = "subject"}}{ groups samples by the target of
  \code{sample_belongs_to_subject}. All samples linked to the same
  \code{Subject} receive the same \code{group_id}.}

  \item{\code{mode = "batch"}}{ groups samples by the target of
  \code{sample_processed_in_batch}. Samples with no batch assignment are
  retained as singleton unlinked groups and recorded in metadata warnings.}

  \item{\code{mode = "study"}}{ groups samples by the target of
  \code{sample_from_study}.}

  \item{\code{mode = "time"}}{ groups samples by the target of
  \code{sample_collected_at_timepoint}. When \code{Timepoint} nodes have
  \code{time_index} metadata, that value is used to derive \code{order_rank}.
  If \code{time_index} is unavailable, the function attempts to derive ordering
  from \code{timepoint_precedes} edges over the timepoint subgraph.}

  \item{\code{mode = "composite"}, \code{strategy = "strict"}}{ projects the
  selected dependency relations onto a sample graph and assigns one
  \code{group_id} per connected component. This is the transitive-closure
  interpretation of composite dependency grouping.}

  \item{\code{mode = "composite"}, \code{strategy = "rule_based"}}{ evaluates
  dependency assignments in deterministic priority order and groups each sample
  by the highest-priority available dependency source. Lower-priority available
  dependencies are retained in the explanation field.}
}

The returned \code{split_constraint$sample_map} always contains
\code{sample_id}, \code{sample_node_id}, \code{group_id},
\code{constraint_type}, \code{group_label}, and \code{explanation}. Time-aware
constraints also include \code{time_index}, \code{timepoint_id}, and
\code{order_rank} when available.
}
