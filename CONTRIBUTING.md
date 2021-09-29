# Contributing guidelines

> **These guidelines are under construction!**

Reference Code of Conduct

Overview of ways to contribute

Links to GitHub issues and NDCN contact email

## Contributing changes via GitHub

Links to GitHub tutorials for these tasks

- fork this repository
- make changes and commit
- submit a PR with a meaningful message explaining what your change accomplishes and why

## Website, documentation, and tutorials

The [`omicser` website](https://ergonyc.github.io/omicser/)
is built using [`pkgdown`](https://pkgdown.r-lib.org/).

When applicable, please edit the `Rmd` version of a file (rather than the `md`) version.
If you would like to recommend changes to the tutorials listed under the "Articles"
drop-down menu,
you can find the source files in the `vignettes` directory.

The site will need to be rebuilt as follows for these changes to be visible on the website:

```r
library(pkgdown)
pkgdown::build_site()
```

All changes from this command will need to be committed to the repository for the website to be updated.

If you are proposing changes to the front page of the website,
you should edit the top-level `README.Rmd` file in the repository.
The file needs to be knit (to create `README.md`)
and then the site rebuilt for these changes to be published to the website.

## Code fixes and enhancements to the web app

TBA
