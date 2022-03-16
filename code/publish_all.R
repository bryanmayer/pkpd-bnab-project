# publishes full pipeline
wflow_publish(
  c(
    "analysis/_site.yml",
    "analysis/index.Rmd",
    "analysis/opt-sims-background.Rmd",
    "analysis/sensitivity-analysis.Rmd",
    "analysis/bispecific.Rmd",
    "analysis/empirical-case-study.Rmd",
    "analysis/hill-slope-meta.Rmd",
    "analysis/titer-analysis.Rmd"
  )
)
