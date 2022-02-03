# Contributing guidelines

We're glad you're interested in contributing to this project!
Before getting started,
please review our [Code of Conduct](CODE_OF_CONDUCT.md).
We expect all interactions related to this project follow these guidelines.

These contribution guidelines assume you have a [GitHub account](https://github.com/)
and are able to contribute via the project's [GitHub repository](https://github.com/ndcn/omicser).
If you are not familiar with GitHub,
you may find it useful to read through GitHub's [Quickstart guide](https://docs.github.com/en/get-started/quickstart/hello-world).
You can contribute in two main ways via GitHub:

- [Issues](#issues): sharing bug reports and ideas for new features; requires minimal coding and/or GitHub skills
- [Pull requests](#pull-requests): providing updates to the code and/or documentation in the repository; experience with GitHub (and potentially code) suggested

Additional guidance on these two contribution methods are included below.

## Issues

[Issues in GitHub](https://docs.github.com/en/issues/tracking-your-work-with-issues/about-issues) are a great way to communicate with us about this project,
whether it is filing a bug report (about something that doens't quite work for you in the browser)
or suggesting a new feature (adding new functionality to the browser).
Browse the [existing issues](https://github.com/ndcn/omicser/issues) to see whether there is a current conversation related to your idea or problem.

You are welcome to add feedback to existing issues,
including comments and reactions to posts in the issue discussion.
If no conversation about your idea currently exists,
please [submit a new issue](https://github.com/ndcn/omicser/issues/new).
If you are filing an issue about a problem with rendering in the browser,
we appreciate you sharing a screen shot by [attaching a file to the conversation](https://docs.github.com/en/github/writing-on-github/working-with-advanced-formatting/attaching-files)!

## Pull requests

[Pull requests (PRs)]((https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests)) are the GitHub mechanism that allows you to suggest changes to someone else's repository (in this case, our project!).
This section includes three types of information:

- [PR workflow](#pr-workflow): the process to create and submit a PR
- [Documentation PRs](#documentation-prs): specific information about PRs involving documentation (including the project's website and tutorials)
- [Code PRs](#code-prs): code style and format guidelines for PRs involving `shiny` and/or other scripts

### PR workflow

If your PR is correcting a small typo in documentation,
or if you are learning how PRs work,
the most straightforward approach is to use the GitHub web interface to
make changes and submit the PR.
The general process is outlined below,
but is also described in more detail in the [GitHub Quickstart section on GitHub flow](https://docs.github.com/en/get-started/quickstart/github-flow):

1. [Fork](https://github.com/ndcn/omicser/fork) the project
2. Create a new [branch](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-and-deleting-branches-within-your-repository#creating-a-branch).
3. Alter the content on the page using the [editor in the web interface](https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files).
4. [Submit a pull request](https://github.com/ndcn/omicser/compare) to our repository from your fork's branch containing the changes.

A member of the core team will review your PR and inform you how to proceed. (see also our [GOVERNANCE document](https://github.com/ndcn/omicser/GOVERNANCE_draft.md))

If you are making more complex changes involving multiple files and/or code,
it is likely you will need to make your changes on a local version of the repository.
The PR workflow then requires Git to be installed and accessible on your local machine,
such as the [GitHub Desktop app](https://docs.github.com/en/desktop), 
[GitHub CLI](https://docs.github.com/en/github-cli),
or [git command line tools](https://git-scm.com/book/en/v2/Getting-Started-The-Command-Line).  
Then the workflow becomes:

1. [Fork](https://github.com/ndcn/omicser/fork) the project
2. Create a new [branch](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches).
3. [Clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) the project to create a local copy on your own computer.
4. Alter the code, check to ensure the changes still allow the browser to launch, and commit your changes.
5. [Push](https://docs.github.com/en/get-started/using-git/pushing-commits-to-a-remote-repository) your changes to the remote repository.
6. [Submit a pull request](https://github.com/ndcn/omicser/compare) to our repository from your fork's branch containing the changes.

The specific steps may differ slightly depending on the method with which you interact with Git on your own computer.

### Documentation PRs

The [`omicser` website](https://ndcn.github.io/omicser/)
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

### Code PRs

> **These guidelines are under construction!**
