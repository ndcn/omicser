
repo_zip_url <- "https://github.com/ndcn/omicser/archive/refs/heads/main.zip"
temp <- tempfile()
download.file(repo_zip_url, temp)
allfiles <- unzip(zipfile=temp,list=TRUE)$Name
# extract only the files in examples
unzip(zipfile = temp,
      files = allfiles[startsWith(allfiles,'omicser-main/examples/')],
      exdir="examples",junkpaths = TRUE)
