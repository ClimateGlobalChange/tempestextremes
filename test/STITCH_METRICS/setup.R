#Run this before using stitch_metric_framework.R


req_libs<-c("abind","akima","argparse",
            "ggplot2","gtable","grid",
            "knitr","markdown","ncdf4",
            "ncdf4.helpers","PCICt",
            "reshape2","RNetCDF","rmarkdown")
for (r in req_libs){
  install.packages(r, repos='http://cran.us.r-project.org')
}
