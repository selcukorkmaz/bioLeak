\name{create_nodes}
\alias{create_nodes}
\alias{create_edges}
\title{Create Canonical Node and Edge Tables}
\usage{
create_nodes(data, type, id_col, label_col = NULL, attr_cols = NULL, prefix = TRUE, dedupe = TRUE)

create_edges(data, from_col, to_col, from_type, to_type, relation, attr_cols = NULL,
  allow_missing = FALSE, dedupe = TRUE)
}
\arguments{
  \item{data}{A \code{data.frame} containing entity or relationship columns.}
  \item{type, from_type, to_type}{Supported node types such as \code{"Sample"} or \code{"Subject"}.}
  \item{id_col}{Column containing the source identifier for the node type.}
  \item{label_col}{Optional column used for node labels.}
  \item{attr_cols}{Optional columns stored in the \code{attrs} list-column.}
  \item{prefix}{If \code{TRUE}, prepend typed prefixes such as \code{sample:} to node identifiers.}
  \item{dedupe}{If \code{TRUE}, collapse duplicate identifiers or duplicate edges.}
  \item{from_col, to_col}{Source and target identifier columns for edge creation.}
  \item{relation}{Canonical edge type.}
  \item{allow_missing}{If \code{TRUE}, drop rows with missing edge endpoints instead of erroring.}
}
\value{
For \code{create_nodes()}, a \code{graph_node_set}. For \code{create_edges()},
a \code{graph_edge_set}.
}
\description{
Build canonical node and edge tables from ordinary metadata frames.
}
\examples{
meta <- data.frame(
  sample_id = c("S1", "S2"),
  subject_id = c("P1", "P2")
)

samples <- create_nodes(meta, type = "Sample", id_col = "sample_id")
edges <- create_edges(
  meta,
  from_col = "sample_id",
  to_col = "subject_id",
  from_type = "Sample",
  to_type = "Subject",
  relation = "sample_belongs_to_subject"
)
}

