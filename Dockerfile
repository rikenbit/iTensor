# Base Image
FROM bioconductor/bioconductor_docker:RELEASE_3_17

# Install R Packages
RUN R -e "BiocManager::install('mixOmics', ask=FALSE); devtools::install_github('rikenbit/iTensor', \
    upgrade='always', force=TRUE, INSTALL_opts = '--install-tests');\
    tools::testInstalledPackage('iTensor')"