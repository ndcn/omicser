---
title: "Installation"
author: "andy henrie"
date: "9/14/2021"
output:
  html_document: 
    keep_md: yes
    toc: true
  md_document:
    variant: markdown_github
---



### NDCN Browser Installation

Now that we've set up our enviroment, installing the NDCN browser with `devtools::install_github` is easy.  

> NOTE: we will soon add releases which will be the optimal way to install. Stay Tuned.


```r
devtools::install_github("ergonyc/omicser")

```

Thats it!  Now comes the crucial (and fun!) part: [curating](03_data_curation.md) your data.

