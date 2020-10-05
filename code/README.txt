Для подготовки запуска R кода нужно:

Установить R, Rtools, Rstudio

Для полной уверенности в Rtools прописать в консоли Rstudio:

writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

и перезапустить

Прописать следующие команды в консоли в Rstudio:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("preprocessCore")

install.packages(c("tseries", "outliers", "CLSOCP", "nortest", "quadprog", "zoo", "fields", "fpc", "flexclust", "plyr", "fastcluster", "skmeans", "LICORS"))
install.packages("C:/Users/nogay/Desktop/gene expression/code/CONORData_1.0.2.tar.gz", repos = NULL, type="source")
install.packages("C:/Users/nogay/Desktop/gene expression/code/HARMONY.tar.gz", repos = NULL, type="source")