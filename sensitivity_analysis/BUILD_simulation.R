### BUILD_simulation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2022 (09:44) 
## Version: 
## Last-Updated: dec  1 2022 (09:23) 
##           By: Brice Ozenne
##     Update #: 44
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## cd /projects/biostat01/people/hpl802/pipeline/
## * Path
if(system("whoami",intern=TRUE) == c("unicph\\hpl802")){
    path <- "x:/pipeline/"
}else if(system("whoami",intern=TRUE) == c("hpl802")){
    path <- ""
}
path.results <- file.path(path,"Results") ## path.results <- "Results"
path.results1 <- file.path(path.results,"simulation-scenario1")
path.results2 <- file.path(path.results,"simulation-scenario2")
path.results3 <- file.path(path.results,"simulation-scenario3")

## * Dependencies
library(data.table)
library(pbapply)

loadRes <- function(path, tempo.file = FALSE, type = NULL,
                    export.attribute = NULL, trace = TRUE){
    all.files <- list.files(path)
    file.tempo <- grep("(tempo)",all.files,value = TRUE)
    file.final <- setdiff(all.files, file.tempo)

    if(tempo.file){
        file.read <- file.tempo
    }else{
        file.read <- file.final
    }
    if(!is.null(type)){
        file.read <- grep(pattern=type,x=file.read,value=TRUE)
    }

    n.file <- length(file.read)

    myApply <- switch(as.character(as.logical(trace)),
                      "TRUE" = pbapply::pblapply,
                      "FALSE" = lapply)

    ls.out <- do.call(myApply, args = list(X = 1:n.file, FUN = function(iFile){
        iRead <- try(readRDS(file = file.path(path,file.read[iFile])))
        if(inherits(iRead,"try-error")){
            return(NULL)
        }else{
            iOut <- cbind(data.table::as.data.table(iRead),
                          file = file.read[iFile])
        return(iOut)
        }
    }))
    out <- do.call(rbind, ls.out)
    return(out)
}

## * Load results
dt.sim1 <- loadRes(path.results1, tempo.file = TRUE)
dt.sim2 <- loadRes(path.results2, tempo.file = TRUE)
dt.sim3 <- loadRes(path.results3, tempo.file = TRUE)
## dt.sim1[n == 1000 & type=="average" & beta == 0, .N, by = c("type","n","beta","file")]

dt.sim <- rbind(cbind(scenario = "scenario 1", dt.sim1),
                cbind(scenario = "scenario 2", dt.sim2),
                cbind(scenario = "scenario 3", dt.sim3))

## * Check
## unique(dt.sim$n)
## dt.sim1[n==1000 & beta == 0 & type %in% c("average","gls"), .(sd(estimate),mean(se)), by = "type"]
## dt.sim2[n==1000 & beta == 0 & type %in% c("average","gls"), .(sd(estimate),mean(se)), by = "type"]
## dt.sim3[n==1000 & beta == 0 & type %in% c("average","gls"), .(sd(estimate),mean(se)), by = "type"]

## * Process results
dt.sim[, n.char := factor(n/2, unique(n/2))]
dt.sim[, beta.char := ifelse(beta==0,paste0("Null hypothesis (\u03B2=",beta,")"),paste0("Alternative hypothesis (\u03B2=",beta,")"))]
dt.sim[, beta.char := factor(beta.char,levels=unique(beta.char))]
dt.sim[, type.char := factor(type,levels=c("average","fixse","gls","gls1","robust","proportion"),c("pool (average)","pool (se)","pool (gls)", "pool (constrained gls)","pool (robust gls)","proportion"))]
dt.sim[, target := factor(type,levels=c("average","fixse","gls","gls1","robust","proportion"),c("pool","pool","pool","pool","pool","proportion"))]
dt.sim[, GS := beta]
dt.sim[type=="proportion", GS := as.numeric(beta>0)]
dtS.sim <- dt.sim[,.(rep = .N, average = mean(estimate), sd = sd(estimate), average.se = mean(se), bias = mean(estimate)-GS, power = mean(p.value<=0.05)),
                    by = c("scenario","beta","n","type","n.char","beta.char","type.char","target","GS")]

## * Export results
saveRDS(dt.sim, file.path(path.results,"simulation-scenario.rds"))
saveRDS(dtS.sim, file.path(path.results,"simulation-summary-scenario.rds"))

##----------------------------------------------------------------------
### BUILD_simulation.R ends here
