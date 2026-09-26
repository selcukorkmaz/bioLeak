# Capture session provenance information

Internal helper that records R version, loaded packages, platform,
timestamp, and optionally git SHA and hardware details.

## Usage

``` r
.bio_capture_provenance(capture_git = TRUE, capture_hardware = TRUE)
```

## Arguments

- capture_git:

  Logical; if TRUE, attempts to capture the current git short SHA via
  `system2("git", ...)`. Default FALSE.

- capture_hardware:

  Logical; if TRUE, captures basic hardware info from
  [`Sys.info()`](https://rdrr.io/r/base/Sys.info.html). Default TRUE.

## Value

A list with elements `r_version`, `packages` (data.frame), `platform`,
`timestamp`, `git_sha`, and `hardware`.
