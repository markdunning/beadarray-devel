library(GEOmetadb)
getSQLiteFile()


con <- dbConnect(SQLite(),'GEOmetadb.sqlite')

sFiles <- geoConvert("GPL6947")$sMatrix
geo_tables <- dbListTables(con)

head(algse)

which.max(lapply(unlist(split(algse[,2], algse[,1])),length))

gsegpl <- dbGetQuery(con, "select * from gse_gpl")
gsegsm <- dbGetQuery(con, "select * from gse_gsm")

illdat <- gsegpl[gsegpl[,2] == "GPL6947",]
sFile <- sapply(illdat[,1], function(x) sFiles[grep(x, sFiles[,2]),2])
which.max(lapply(split(illdat[,2],illdat[,1]),length))
