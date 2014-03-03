test_that("Vector output format", {
            
            # MM vector output
            file <- system.file("extdata", "test1.csv", package="PhenStat")
            test <- PhenList(dataset=read.csv(file),
                    testGenotype="Sparc/Sparc")
            result <- testDataset(test,
                    depVariable="Lean.Mass")
            vector_results <- vectorOutput(result)    
            expect_that(length(vector_results), equals(32))    
            expect_that(is.null(names(vector_results)), is_false())       
            expect_that(length(names(vector_results)), equals(32))    
            expect_that(all(lapply(vector_results,is.character)==TRUE),is_true())
            
            # FE vector output
            file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
            test2 <- PhenList(dataset=read.csv(file),
                    testGenotype="Aff3/Aff3")
            result2 <- testDataset(test2,
                    depVariable="Thoracic.Processes",
                    method="FE")  
            vector_results2 <- vectorOutput(result2)
            expect_that(length(vector_results2), equals(32))    
            expect_that(is.null(names(vector_results2)), is_false())       
            expect_that(length(names(vector_results2)), equals(32)) 
            expect_that(all(lapply(vector_results2,is.character)==TRUE),is_true())
            
            # RR vector output
            # is coming    
                 
}
)