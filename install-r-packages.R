packages_file <- "r-packages.txt"

if (!file.exists(packages_file)) {
  stop("Cannot find ", packages_file, ". Run this script from the eu-aiv-analysis project folder.")
}

packages <- readLines(packages_file, warn = FALSE)
packages <- trimws(packages)
packages <- packages[nzchar(packages)]
packages <- unique(packages)

installed <- rownames(installed.packages())
missing <- setdiff(packages, installed)

if (length(missing) == 0) {
  message("All R packages in ", packages_file, " are already installed.")
} else {
  message("Installing missing R packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}
