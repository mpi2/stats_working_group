Example of usage:
> getIMPCTable("./IMPCData.csv","WTSI","MGP_001","IMPC_CBC_001")

Results are saved in IMPCData.csv file.

Take the first row and the last column with dataset function call:
> IMPC_dataset1 <- getIMPCDataset('WTSI','MGP_001','IMPC_CBC_001','IMPC_CBC_003_001','MGI:4431644')

Now dataset is ready to use with PhenStat:

> testIMPC1 <- PhenList(dataset=IMPC_dataset1,
        testGenotype="MDTZ",refGenotype="+/+",
        dataset.colname.genotype="Colony")

> resultIMPC1_RR <- testDataset(testIMPC1,
        depVariable="Value",
        method="RR")
> summaryOutput(resultIMPC1_RR)

> resultIMPC1_MM <- testDataset(testIMPC1,
        depVariable="Value",
        method="MM")
> summaryOutput(resultIMPC1_MM)