# vcr configuration — runs once when testthat loads helpers
library(vcr)

vcr_dir <- here::here("tests", "fixtures")

vcr_configure(
  dir             = vcr_dir,
  record          = "none",       # "none" in CI; change to "new_episodes" to add cassettes
  log             = TRUE,
  log_opts        = list(file = file.path(vcr_dir, "vcr.log")),
  # Allow real requests when no cassette is active (e.g. in interactive use)
  allow_passthrough = !nzchar(Sys.getenv("CI"))
)
