
library(ogttMetrics)

context("data manipulation")

data(omnicBase_samp)
data(omnicbCbG_samp) # omnicCG

test_that("omnicBase_samp properly loads", {
 expect_true(length(experiments(omnicBase_samp))==2)
 expect_equal(dim(colData(omnicBase_samp)),c(40,4))
 expect_equal(dim(experiments(omnicBase_samp)$glucose_base),c(7,40))
})

test_that("subsetByAssay succeeds", {
 lit = subsetByAssay(omnicBase_samp, c("glucose_base", "insulin_base"))
 expect_true(length(experiments(lit))==2)
})

test_that("addAssay succeeds", {
 lit = subsetByAssay(omnicBase_samp, c("glucose_base", "insulin_base"))
 d = dim(ee <- experiments(omnicBase_samp)$glucose_base)
 dum = matrix(rnorm(prod(d)), nrow=d[1])
 rownames(dum) = rownames(ee)
 colnames(dum) = colnames(ee)
 lit = addAssay(lit, dum, "dum")
 expect_true(length(experiments(lit))==3)
})

context("matsuda")

test_that("matsuda120 agrees with web calculator", {
 jj =  ogttMetrics:::matsuda120( c(100, 160, 160, 160, 140), c(5, 10, 10, 10, 5), c(0, 30, 60, 90, 120)) -> jj
 expect_true(round(jj,2) == 12.34)
})
test_that("addMatsuda120 succeeds", {
 lit = subsetByAssay(omnicBase_samp, c("glucose_base", "insulin_base"))
 lit = addMatsuda120(lit, gname="glucose_base", iname="insulin_base")
 expect_true(length(experiments(lit))==3)
 expect_true(abs(experiments(lit)$Mats120[1] - 14.62961) < 1e-4)
})

context("minmod")

test_that("getMinmodSIs succeeds", {
 lit = omnicBase_samp[,1:5]
 bw = colData(lit)[,"bodyweight"]
 bw = matrix(bw,nrow=1)
 rownames(bw) = "bodyweight"
 colnames(bw) = colnames(lit)[[1]]
 lit = addAssay(lit, bw, "bodyweight")
 suppressWarnings({sis = getMinmodSIs(lit,
    gname="glucose_base", iname="insulin_base")}) # some known non-convergence for these data
 expect_true(all(dim(sis)==c(5,5)))
 expect_true(abs(round(sis[1,3],4)-0.0011)<.0001)
})

context("plot_OGTT_fit")

test_that("plot_OGTT_fit does not error", {
 tf = tempfile()
 pdf(tf)
 bw = colData(omnicBase_samp)[,"bodyweight"]
 bw = matrix(bw,nrow=1)
 rownames(bw) = "bodyweight"
 colnames(bw) = colnames(omnicBase_samp)[[1]]
 omnicBase_samp = addAssay(omnicBase_samp, bw, "bodyweight")
 suppressWarnings({mm = minmodByID(omnicBase_samp, "1",
   iname="insulin_base", gname="glucose_base")})
 pl <- plot_OGTT_fit(mm)
 pl
 dev.off()
 expect_true(file.exists(tf))
 expect_true(class(pl)[1] == "gg")
})

context("concatenation of MAEs")

test_that("concat is reliable", {
 tstcon = concat(list(base=omnicBase_samp, CG=omnicCG_samp),
    basenames=c("glucose", "insulin"))
 expect_true(length(experiments(tstcon))==4)
 expect_equal(names(experiments(tstcon)), c("glucose_base", "insulin_base", "glucose_CG", "insulin_CG"))
# expect_equal(experiments(tstcon)$Mats120_base[1:5],
#     experiments(omnicBase_samp)$Mats120[1:5])
})

context("outlier detection")

test_that("mvOutliers succeeds", {
 data(obaSamp)
 mm = mvOutliers(obaSamp)
 expect_true(mm$insulinOutliers == "34")
 expect_true(is.na(mm$glucoseOutliers))
})

context("csv import")

test_that("csvImport obtains correct entity", {
 pref = system.file("csv_example", package="ogttMetrics")
 gpath = paste0(pref, "/glucBase.csv")
 ipath = paste0(pref, "/insBase.csv")
 spath = paste0(pref, "/sampBase.csv")
 democ = csvImport(gpath, ipath, spath)
 expect_true(class(democ)=="ogttCohort")
 expect_true(all(names(rownames(democ))==c("glucose", "insulin")))
})

context("QCplots")

test_that("QCplots does not fail", {
 data(obaSamp)
 tf = tempfile()
 pdf(tf)
 QCplots(obaSamp)
 dev.off()
 expect_true(file.exists(tf))
})

context("addMinmodSIs")

test_that("addMinmodSIs errors appropriately", {
 data(obaSamp)
 expect_error(addMinmodSIs(obaSamp[,1:3], replace=FALSE))
})
