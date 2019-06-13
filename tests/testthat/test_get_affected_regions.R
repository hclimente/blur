net <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)

cones <- minigwas$map
colnames(cones) <- c("chr","snp","cm","pos","allele.1", "allele.2")
cones$selected <- FALSE
cones$c <- abs(rnorm(nrow(cones), 0, 1))

# region 1, spanning genes A and B
cones$selected[4:8] <- TRUE
cones$c[4:8] <- abs(rnorm(5, 10, 1))
# region 2, spanning gene C
cones$selected[16:21] <- TRUE
cones$c[16:21] <- abs(rnorm(6, 10, 1))

cones <- martini:::get_snp_modules(cones, net)

regions <- get_affected_regions(cones, net)

test_that("retrieve the right regions", {

  expected_regions <- read.table(text = "
                       module chr start end genes
                       1 1 40 70 A,B
                       2 2 45 85 C
                       ", header = TRUE, stringsAsFactors = FALSE)

  expect_equal(regions, expected_regions)
})
