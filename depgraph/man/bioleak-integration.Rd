\name{as_bioleak_split_spec}
\alias{as_bioleak_split_spec}
\alias{validate_bioleak_split_spec}
\alias{summarize_leakage_risks}
\alias{bioleak_split_spec}
\alias{bioleak_split_spec_validation}
\alias{leakage_risk_summary}
\title{Translate depgraph Constraints for bioLeak}
\usage{
as_bioleak_split_spec(constraint, graph = NULL)

validate_bioleak_split_spec(x)

summarize_leakage_risks(graph, constraint = NULL, split_spec = NULL,
  validation = NULL)
}
\arguments{
  \item{constraint}{A \code{split_constraint}.}
  \item{graph}{A \code{dependency_graph}.}
  \item{split_spec}{An optional \code{bioleak_split_spec}.}
  \item{validation}{An optional \code{depgraph_validation_report}.}
  \item{x}{A \code{bioleak_split_spec}.}
}
\value{
\code{as_bioleak_split_spec()} returns a \code{bioleak_split_spec}.
\code{validate_bioleak_split_spec()} returns a
\code{bioleak_split_spec_validation}. \code{summarize_leakage_risks()}
returns a \code{leakage_risk_summary}.
}
\description{
Translate graph-derived split constraints into a stable, inspectable structure
that downstream \code{bioLeak} workflows can consume, perform preflight
structural checks on that translation, and summarize structural leakage risks.
}
\details{
The translation layer always produces canonical sample-level columns including
\code{sample_id}, \code{sample_node_id}, \code{group_id}, and
\code{primary_group}. When available, it also carries \code{batch_group},
\code{study_group}, \code{timepoint_id}, \code{time_index}, and
\code{order_rank}. Missing but relevant fields are retained as \code{NA}
columns rather than omitted.

The split-spec validator checks:
\itemize{
  \item missing required columns
  \item missing or duplicated sample identifiers
  \item missing grouping assignments
  \item singleton-only grouping structures
  \item missing ordering when ordering is required
  \item invalid or empty block variables
}

\code{summarize_leakage_risks()} reuses \code{validate_graph()} and
\code{split_constraint} metadata rather than duplicating \code{bioLeak}
analytics.
}
