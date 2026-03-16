\name{ingest_metadata}
\alias{ingest_metadata}
\title{Standardize Sample Metadata}
\usage{
ingest_metadata(data, col_map = NULL, dataset_name = NULL, strict = TRUE)
}
\arguments{
  \item{data}{A sample-level \code{data.frame}.}
  \item{col_map}{Optional named character vector mapping canonical names to user-provided columns.}
  \item{dataset_name}{Optional dataset label stored as an attribute on the returned table.}
  \item{strict}{If \code{TRUE}, error when required columns are missing.}
}
\value{
A standardized \code{data.frame} with canonical identifier columns coerced to character.
}
\description{
Normalize user-provided metadata into the canonical column contract used by
\pkg{depgraph}.
}
\examples{
meta <- ingest_metadata(
  data.frame(sample_id = c("S1", "S2"), subject_id = c("P1", "P2"))
)
}

